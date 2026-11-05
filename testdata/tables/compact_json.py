"""compact_json.py: copy JSON from stdin to stdout without the whitespace
between tokens (what Go's encoding/json.Compact does). Every token is kept
byte for byte (numbers are not re-formatted), so the compact golden holds
exactly what the cfitsio harness printed. Plain Python, no libraries.

  ../cfref/cfref dump x.fits | python3 compact_json.py > golden/x.fits.json
"""

import sys


def compact(src: bytes) -> bytes:
    out = bytearray()
    in_str = esc = False
    for b in src:
        if in_str:
            out.append(b)
            if esc:
                esc = False
            elif b == 0x5C:  # backslash
                esc = True
            elif b == 0x22:  # closing quote
                in_str = False
        elif b == 0x22:
            in_str = True
            out.append(b)
        elif b not in b" \t\r\n":
            out.append(b)
    return bytes(out)


if __name__ == "__main__":
    sys.stdout.buffer.write(compact(sys.stdin.buffer.read()))
