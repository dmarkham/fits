"""cmp_golden.py COMMITTED REGENERATED: exit 1 unless both directory trees
hold the same files with the same bytes. A <name>.gz file (a gzipped
golden, see put_golden in gen.py) is compared after decompression: its
content is the harness output; the compressed bytes depend on the zlib
build."""
import gzip
import os
import sys


def files(root):
    out = set()
    for d, _, fs in os.walk(root):
        for f in fs:
            out.add(os.path.relpath(os.path.join(d, f), root))
    return out


def content(path):
    b = open(path, "rb").read()
    return gzip.decompress(b) if path.endswith(".gz") else b


def main():
    a, b = sys.argv[1], sys.argv[2]
    fa, fb = files(a), files(b)
    bad = 0
    for f in sorted(fa ^ fb):
        print(f"only in {a if f in fa else b}: {f}")
        bad = 1
    for f in sorted(fa & fb):
        if content(os.path.join(a, f)) != content(os.path.join(b, f)):
            print(f"differs: {os.path.join(a, f)}")
            bad = 1
    sys.exit(bad)


if __name__ == "__main__":
    main()
