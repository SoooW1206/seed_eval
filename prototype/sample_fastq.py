#!/usr/bin/env python3
"""Exactly N reads drawn uniformly at random, without replacement, from a FASTQ on stdin, written to stdout
in their original order.

Usage: pigz -dc IN.fastq.gz | python3 -I sample_fastq.py --total T --n N --seed S --names NAMES.txt | pigz > OUT

--total is the number of reads in the input (known beforehand); the record numbers to keep are
random.Random(S).sample(range(T), N), so the sample depends only on T, N and S. Exits with an error if the
input does not hold exactly T four-line records. NAMES.txt gets the kept read names (first word of the
header, without '@'), one per line. Prints only counts to stderr, never sequences.
"""
import argparse
import random
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--total", type=int, required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--names", required=True)
    args = parser.parse_args()
    if not 0 < args.n <= args.total:
        sys.exit("need 0 < --n <= --total")
    keep = set(random.Random(args.seed).sample(range(args.total), args.n))
    source, out = sys.stdin.buffer, sys.stdout.buffer
    index = written = 0
    with open(args.names, "w") as names:
        while True:
            header = source.readline()
            if not header:
                break
            record = [header, source.readline(), source.readline(), source.readline()]
            if not header.startswith(b"@") or not record[2].startswith(b"+") or not record[3]:
                sys.exit(f"record {index}: not a four-line FASTQ record")
            if index in keep:
                out.write(b"".join(record))
                names.write(header[1:].split()[0].decode() + "\n")
                written += 1
            index += 1
    if index != args.total:
        sys.exit(f"input has {index} reads, --total says {args.total}")
    if written != args.n:
        sys.exit(f"wrote {written} reads, expected {args.n}")
    print(f"reads in input: {index}; kept: {written}; seed: {args.seed}", file=sys.stderr)


if __name__ == "__main__":
    main()
