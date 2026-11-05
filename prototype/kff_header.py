#!/usr/bin/env python3
"""Print k, max, min and data_size from the first global-variable section of a KFF file.

Usage: python3 -I kff_header.py FILE.kff

Reads only the header and the first 'v' section (skipping a leading 'i' index section);
never reads k-mer data. Layout as in vg's deps/kff-cpp-api/kff_io.cpp (all integers big-endian):
header = "KFF", major, minor, encoding, uniqueness, canonicity (1 byte each), metadata size (4 bytes),
metadata; 'v' section = 'v', variable count (8 bytes), then per variable a NUL-terminated name and
an 8-byte value; 'i' section = 'i', entry count (8 bytes), entries of 1 + 8 bytes, next index (8 bytes).
Prints numbers only. A variable missing from the first 'v' section is printed as NA.
"""
import struct
import sys

VARIABLES = ["k", "max", "min", "data_size"]
MAX_HEADER_BYTES = 1 << 20  # header + first variable section are a few hundred bytes in practice


class Reader:
    def __init__(self, f):
        self.f = f
        self.used = 0

    def read(self, n):
        self.used += n
        if self.used > MAX_HEADER_BYTES:
            sys.exit(f"no variable section within the first {MAX_HEADER_BYTES} bytes; stopping")
        data = self.f.read(n)
        if len(data) != n:
            sys.exit("file ends inside the header")
        return data

    def u64(self):
        return struct.unpack(">Q", self.read(8))[0]


def first_variables(path):
    with open(path, "rb") as f:
        r = Reader(f)
        if r.read(3) != b"KFF":
            sys.exit("not a KFF file (no KFF signature)")
        r.read(5)  # major, minor, encoding, uniqueness, canonicity
        r.read(struct.unpack(">I", r.read(4))[0])  # metadata
        while True:
            section = r.read(1)
            if section == b"i":
                r.read(r.u64() * 9 + 8)
            elif section == b"v":
                variables = {}
                for _ in range(r.u64()):
                    name = bytearray()
                    while (c := r.read(1)) != b"\0":
                        name += c
                    variables[name.decode("ascii", "replace")] = r.u64()
                return variables
            else:
                sys.exit("first data section comes before any variable section")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    variables = first_variables(sys.argv[1])
    for name in VARIABLES:
        print(f"{name}\t{variables.get(name, 'NA')}")


if __name__ == "__main__":
    main()
