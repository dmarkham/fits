#!/usr/bin/env bash
# Run the test suite on a Mac from here. The repo is edited on this
# machine; each command copies the working tree to the Mac (respecting
# .gitignore, so no Linux harness binaries or venv go across) and runs there
# over SSH.
#
#   tools/mac.sh check     # macOS, CPU, Go, Homebrew and the reference tools
#   tools/mac.sh setup     # one-time: Homebrew go, cfitsio, wcslib, wcstools, pyenv; venv in ~/.venvs/fits
#   tools/mac.sh test      # tier 1: go vet, then go test ./... with PATH limited to Go's bin dir
#   tools/mac.sh interop   # tier 2: build the C harnesses, then go test -tags interop ./...
#   tools/mac.sh run CMD   # run CMD in the synced tree (e.g. 'go test -run TestLock -v .')
#
# MAC_HOST is an ssh host (alias in ~/.ssh/config), default "mac".
# REMOTE_DIR is the copy on the Mac, default src/fits.
# Extra arguments to test/interop are passed to go test (e.g. -run TestEdit -v).
set -euo pipefail

cd "$(dirname "$0")/.."
MAC_HOST=${MAC_HOST:-mac}
REMOTE_DIR=${REMOTE_DIR:-src/fits}
cmd=${1:-test}
shift || true

sync() {
  ssh "$MAC_HOST" "mkdir -p '$REMOTE_DIR'"
  # Commits first, as a git bundle, so the Mac copy is a real checkout with
  # history (some tier 2 tests read git). Only the commits the Mac lacks
  # are sent once it has a history.
  local have bundle
  have=$(ssh "$MAC_HOST" "git -C '$REMOTE_DIR' rev-parse -q --verify HEAD 2>/dev/null" || true)
  bundle=$(mktemp "${TMPDIR:-.}/fits-mac.XXXXXX.bundle")
  if [[ -n "$have" ]] && git merge-base --is-ancestor "$have" HEAD 2>/dev/null; then
    [[ "$have" == "$(git rev-parse HEAD)" ]] || git bundle create -q "$bundle" "$have..HEAD"
  else
    git bundle create -q "$bundle" HEAD
  fi
  if [[ -s "$bundle" ]]; then
    scp -q "$bundle" "$MAC_HOST:$REMOTE_DIR.bundle"
    ssh "$MAC_HOST" "cd '$REMOTE_DIR' && { [ -d .git ] || git init -q; } && git fetch -q '../$(basename "$REMOTE_DIR").bundle' HEAD && git checkout -q -f --detach FETCH_HEAD && rm -f '../$(basename "$REMOTE_DIR").bundle'"
  fi
  rm -f "$bundle"
  # Then the working tree (uncommitted edits included), retrying because the
  # Mac's connection sometimes drops mid-transfer. Gitignored files are not
  # sent, and --delete may remove them on the Mac, so the Mac's venv lives
  # outside the tree ($VENV) and `interop` rebuilds the C harnesses.
  local i
  for i in 1 2 3; do
    rsync -az --partial --delete --exclude .git --exclude /venv/ --filter=':- .gitignore' \
      ./ "$MAC_HOST:$REMOTE_DIR/" && return
    echo "rsync failed; retrying" >&2
    sleep 3
  done
  return 1
}

# Runs the remote half; everything below the marker executes on the Mac.
remote() {
  # caffeinate keeps the laptop from idle-sleeping mid-command.
  ssh -o ServerAliveInterval=15 "$MAC_HOST" "caffeinate -i bash -s -- $(printf '%q ' "$@")" <<'REMOTE'
set -euo pipefail
cmd=$1 dir=$2
shift 2
[[ -x /opt/homebrew/bin/brew ]] && eval "$(/opt/homebrew/bin/brew shellenv)"
[[ -x /usr/local/bin/brew ]] && eval "$(/usr/local/bin/brew shellenv)"
cd "$dir" 2>/dev/null || true
VENV=$HOME/.venvs/fits
if command -v brew >/dev/null; then
  # The harness Makefiles use plain -lcfitsio/-lwcs and <wcslib/...>.
  export CPATH="$(brew --prefix)/include${CPATH:+:$CPATH}"
  export LIBRARY_PATH="$(brew --prefix)/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
fi

case "$cmd" in
  check)
    sw_vers
    echo "arch: $(uname -m)"
    for t in brew go cc python3 fpack funpack fitsverify xy2sky pkg-config; do
      printf '%-11s %s\n' "$t" "$(command -v $t || echo MISSING)"
    done
    command -v go >/dev/null && go version
    command -v brew >/dev/null && for f in cfitsio wcslib; do
      printf '%-11s %s\n' "$f" "$(brew list --versions $f || echo 'not installed')"
    done
    [[ -x $VENV/bin/python ]] && $VENV/bin/python -c 'import numpy, astropy; print("venv: numpy", numpy.__version__, "astropy", astropy.__version__)' \
      || echo "venv: missing (tools/mac.sh setup)"
    ;;
  setup)
    command -v brew >/dev/null || { echo "install Homebrew first: https://brew.sh" >&2; exit 1; }
    # cfitsio brings fpack, funpack and fitsverify; wcstools brings xy2sky/sky2xy;
    # pkgconf + glib are for stretch/testdata/sref (Siril C sources).
    # openssl@3, readline, xz: what pyenv's python-build needs on macOS.
    for f in go cfitsio wcslib wcstools pkgconf glib pyenv openssl@3 readline xz; do
      brew list --versions $f >/dev/null || brew install $f
    done
    # Python as on the Linux side: pyenv 3.11.11 and a local venv.
    # A build without _ssl is broken (no pip from PyPI); rebuild it.
    py="$(pyenv root)/versions/3.11.11/bin/python"
    if [[ -x "$py" ]] && ! "$py" -c 'import ssl' 2>/dev/null; then pyenv uninstall -f 3.11.11; fi
    (unset CPATH LIBRARY_PATH; pyenv install -s 3.11.11)
    [[ -x $VENV/bin/python ]] && $VENV/bin/python -c 'import ssl' 2>/dev/null || { rm -rf "$VENV"; mkdir -p "$(dirname "$VENV")"; "$py" -m venv "$VENV"; }
    $VENV/bin/pip install -q --upgrade pip
    $VENV/bin/pip install -q numpy==2.4.6 astropy==8.0.2 opencv-python-headless==5.0.0.93
    ;;
  test)
    # Tier 1 needs nothing but Go: run it with every other program out of
    # reach (PATH = Go's bin directory only, no FITS_* variables).
    go vet ./...
    gobin=$(dirname "$(command -v go)")
    env -i HOME="$HOME" PATH="$gobin" GOCACHE="$(go env GOCACHE)" GOPATH="$(go env GOPATH)" \
      GOMODCACHE="$(go env GOMODCACHE)" GOFLAGS= GOTOOLCHAIN=local TMPDIR="${TMPDIR:-/tmp}" \
      go test -count=1 "$@" ./...
    ;;
  interop)
    # Tier 2: the reference tools live. Every missing tool fails its test.
    for mk in $(find . -name Makefile -not -path './venv/*' | sort); do
      d=$(dirname "$mk")
      # These two are used only to regenerate goldens (healpix_cxx and the
      # Siril/OpenCV oracle), which is byte-exact on amd64 only
      # (testutil.RequireGoldenPlatform); their tests skip here.
      case "$d" in ./healpix/testdata/hpref|./stretch/testdata/sref)
        [[ "$(uname -m)" == x86_64 ]] || { echo "skip    $d (golden regeneration only; amd64)"; continue; } ;;
      esac
      make -s -C "$d" >/dev/null 2>&1 && echo "built   $d" || echo "FAILED  $d (its tier 2 tests will fail)"
    done
    go vet -tags interop ./...
    FITS_ASTROPY_PYTHON="$VENV/bin/python" go test -tags interop -count=1 "$@" ./...
    ;;
  run)
    export FITS_ASTROPY_PYTHON="$VENV/bin/python" VENV
    bash -c "$*"
    ;;
  *) echo "unknown command: $cmd" >&2; exit 2 ;;
esac
REMOTE
}

case "$cmd" in
  check) remote "$cmd" "$REMOTE_DIR" ;;
  setup|test|interop|run)
    sync
    remote "$cmd" "$REMOTE_DIR" "$@"
    ;;
  *) echo "usage: tools/mac.sh check|setup|test|interop|run CMD" >&2; exit 2 ;;
esac
