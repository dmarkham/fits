#!/usr/bin/env bash
# Build siril-cli 1.4.0 from the release tarball, for gen_golden.py and
# TestSirilCLIReproduce. Nothing is installed system-wide: everything lives
# under $ROOT (required: a build directory holding the release tarball).
#
# Siril 1.4.0 needs OpenCV 4 and this system has OpenCV 5, so OpenCV 4.13.0
# is the Arch package opencv-4.13.0-9 unpacked (not installed) into
# $ROOT/deps/opencv4, with a pkg-config file whose prefix points there and
# whose Libs list only the modules Siril links (core, imgproc, calib3d and
# what calib3d needs). Optional features Siril does not need here are off:
# ffmpeg (the system ffmpeg is newer than OpenCV 4.13's), libjxl, libheif,
# criterion. libcurl stays on: Siril 1.4.0 does not compile without it
# (core/siril_networking.c uses gboolean without including glib).
#
# Result: $ROOT/siril-cli-1.4.0, a wrapper that sets LD_LIBRARY_PATH to the
# unpacked OpenCV 4 and runs $ROOT/build-cli/src/siril-cli.
#
# Built on 2026-11-01 with gcc 16.2.1, glibc 2.44, meson 1.12.1, ninja 1.13.2.
set -euo pipefail

ROOT=${ROOT:?set ROOT to a build directory holding the Siril 1.4.0 tree}
SIRIL_TGZ_SHA=439def7c40ad783afb82e87f3c656d85b449c701f67ab0a9b97ca372ea2d73c9
OPENCV_PKG=opencv-4.13.0-9-x86_64.pkg.tar.zst
OPENCV_SHA=69ce38e8db804ea58f187996e4771b9a6245c762cd22c047a0ee7cb93c107635

cd "$ROOT"
[ -f siril-1.4.0.tar.bz2 ] || curl -sSfLO https://free-astro.org/download/siril-1.4.0.tar.bz2
echo "$SIRIL_TGZ_SHA  siril-1.4.0.tar.bz2" | sha256sum -c -
[ -d src/siril-1.4.0 ] || { mkdir -p src && tar -C src -xjf siril-1.4.0.tar.bz2; }

mkdir -p deps/pc
[ -f "deps/$OPENCV_PKG" ] || curl -sSfL -o "deps/$OPENCV_PKG" "https://archive.archlinux.org/packages/o/opencv/$OPENCV_PKG"
echo "$OPENCV_SHA  deps/$OPENCV_PKG" | sha256sum -c -
[ -d deps/opencv4/usr/lib ] || { mkdir -p deps/opencv4 && tar -C deps/opencv4 -xf "deps/$OPENCV_PKG"; }
sed -e "s#^prefix=.*#prefix=$ROOT/deps/opencv4/usr#" \
    -e 's#^Libs: .*#Libs: -L${exec_prefix}/lib -lopencv_calib3d -lopencv_features2d -lopencv_flann -lopencv_photo -lopencv_imgproc -lopencv_core#' \
    deps/opencv4/usr/lib/pkgconfig/opencv4.pc > deps/pc/opencv4.pc

if [ ! -d build-cli ]; then
	(cd src/siril-1.4.0 && PKG_CONFIG_PATH="$ROOT/deps/pc" meson setup --buildtype release \
		-Dffmpeg=false -Dlibjxl=false -Dlibheif=false -Dcriterion=false -Dlibcurl=true "$ROOT/build-cli")
fi
ninja -C build-cli src/siril-cli

cat > siril-cli-1.4.0 <<EOF
#!/bin/sh
# siril-cli 1.4.0 (build_siril_cli.sh): release tarball, OpenCV 4.13.0 unpacked locally.
LD_LIBRARY_PATH=$ROOT/deps/opencv4/usr/lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH} exec $ROOT/build-cli/src/siril-cli "\$@"
EOF
chmod +x siril-cli-1.4.0
./siril-cli-1.4.0 --version 2>/dev/null | grep '^siril 1.4.0$'
