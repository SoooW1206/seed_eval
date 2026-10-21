#!/usr/bin/env python3
"""Prototype of minimaxmer selection (Yoshi Suzuki's himapper_rp), steps 1-2 only.

A minimaxmer is chosen like a minimizer, but each k-mer is ranked by the key
    (genome_count ascending, reads_count descending, hash ascending)
instead of by its hash alone: prefer k-mers that are rare in the genome, then
k-mers that are common in the reads. Step 1 selects the minimum of every window
of w consecutive k-mers; step 2 counts how many windows each selected k-mer won.
Yoshi's step 3 (uniformly sparse thinning) is not implemented.

To make a later port to vg/gbwtgraph direct, everything except the key follows
gbwtgraph (deps/gbwtgraph/include/gbwtgraph/minimizer.h, MinimizerIndex::minimizers):
  * 2-bit encoding A=0 C=1 G=2 T=3 (either case), first base in the high bits;
  * hash = wang_hash_64 of the encoding; the canonical orientation of a k-mer is
    the strand with the smaller hash (forward on ties), and the counts are looked
    up with that orientation's code;
  * a k-mer containing any other character is not a candidate;
  * every tied minimum of a window is reported, each read offset at most once.
With constant counts the selection is exactly gbwtgraph's (w, k)-minimizers.
As in Yoshi's code, k-mers with genome_count == 0 are not candidates (skip_absent).

Command line (statistics only, never sequences):
    python3 -I minimaxmer.py --k 31 --w 20 --stats < sequences.txt
reads one sequence per line, counts the k-mers of all of them for reads_count,
uses a constant genome_count of 1, and prints per-read selection counts
(--per-read FILE) and a summary.
"""
import argparse
import collections
import sys

MASK64 = (1 << 64) - 1
PACK = {c: i for i, c in enumerate("ACGT")} | {c: i for i, c in enumerate("acgt")}


def wang_hash_64(key):
    """gbwt::wang_hash_64 on a 64-bit unsigned integer."""
    key = ((key ^ MASK64) + (key << 21)) & MASK64
    key ^= key >> 24
    key = (key + (key << 3) + (key << 8)) & MASK64
    key ^= key >> 14
    key = (key + (key << 2) + (key << 4)) & MASK64
    key ^= key >> 28
    key = (key + (key << 31)) & MASK64
    return key


def encode(kmer):
    """gbwtgraph Key64 code of a k-mer (forward orientation)."""
    code = 0
    for c in kmer:
        code = (code << 2) | PACK[c]
    return code


def reverse_complement_code(code, k):
    rc = 0
    for _ in range(k):
        rc = (rc << 2) | (3 - (code & 3))
        code >>= 2
    return rc


def canonical_kmers(seq, k):
    """For each k-mer start offset: (code, hash, is_reverse) of the canonical orientation, or None
    if the k-mer contains a non-ACGT character. Rolling updates as in gbwtgraph Key64::forward/reverse."""
    low = (1 << (2 * k)) - 1
    forward = reverse = valid = 0
    out = []
    for i, c in enumerate(seq):
        p = PACK.get(c)
        if p is None:
            forward = reverse = valid = 0
        else:
            forward = ((forward << 2) | p) & low
            reverse = ((p ^ 3) << (2 * (k - 1))) | (reverse >> 2)
            valid += 1
        if i >= k - 1:
            if valid >= k:
                fh, rh = wang_hash_64(forward), wang_hash_64(reverse)
                out.append((reverse, rh, True) if rh < fh else (forward, fh, False))
            else:
                out.append(None)
    return out


def select_minimaxmers(seq, k, w, genome_count=None, reads_count=None, skip_absent=True):
    """Steps 1-2. genome_count(code) and reads_count(code) take the canonical code of a k-mer
    (None means a constant count). Returns (selected, wins):
      selected: list of (offset, is_reverse, code), offset = start of the k-mer on the sequence,
                sorted by offset; wins: {offset: number of windows whose minimum it was}."""
    if len(seq) < k + w - 1:
        return [], {}
    kmers = canonical_kmers(seq, k)
    keys = []
    for kmer in kmers:
        if kmer is None:
            keys.append(None)
            continue
        code, h, _ = kmer
        g = genome_count(code) if genome_count else 1
        r = reads_count(code) if reads_count else 0
        keys.append(None if (skip_absent and g == 0) else (g, -r, h))

    # Monotone deque of candidate offsets, as gbwtgraph's CircularBuffer: equal keys stay.
    buffer = collections.deque()
    selected, wins = [], {}
    last_key, last_offset, next_offset = None, -1, 0
    for pos, key in enumerate(keys):
        if buffer and buffer[0] + w <= pos:
            buffer.popleft()
        if key is not None:
            while buffer and keys[buffer[-1]] > key:
                buffer.pop()
            buffer.append(pos)
        if pos < w - 1 or not buffer:
            continue
        front = keys[buffer[0]]
        for i in buffer:
            if keys[i] != front:
                break
            wins[i] = wins.get(i, 0) + 1
        # gbwtgraph's reporting rule.
        if not selected or last_key == front or last_offset < buffer[0]:
            for i in buffer:
                if keys[i] != front:
                    break
                if i >= next_offset:
                    selected.append(i)
                    last_key, last_offset, next_offset = keys[i], i, i + 1
    selected.sort()
    return [(i, kmers[i][2], kmers[i][0]) for i in selected], {i: wins[i] for i in selected}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--w", type=int, default=20)
    parser.add_argument("--per-read", help="write read index, length and number of selected k-mers here")
    args = parser.parse_args()

    sequences = [line.strip() for line in sys.stdin if line.strip()]
    reads = collections.Counter()
    for seq in sequences:
        reads.update(kmer[0] for kmer in canonical_kmers(seq, args.k) if kmer is not None)
    per_read, total_bp, total_selected, total_wins = [], 0, 0, collections.Counter()
    for index, seq in enumerate(sequences):
        selected, wins = select_minimaxmers(seq, args.k, args.w, reads_count=reads.__getitem__)
        per_read.append((index, len(seq), len(selected)))
        total_bp += len(seq)
        total_selected += len(selected)
        total_wins.update(wins.values())
    if args.per_read:
        with open(args.per_read, "w") as f:
            f.write("read_index\tlength\tselected\n")
            f.writelines(f"{i}\t{n}\t{s}\n" for i, n, s in per_read)
    print(f"reads\t{len(sequences)}\ntotal_bp\t{total_bp}\nselected\t{total_selected}")
    print(f"density\t{total_selected / total_bp:.6g}\ntheory_2/(w+1)\t{2 / (args.w + 1):.6g}")
    print(f"distinct_kmers_in_reads\t{len(reads)}\nkmers_seen_once\t{sum(1 for c in reads.values() if c == 1)}")
    print("windows_won\tselected_kmers")
    for won in sorted(total_wins):
        print(f"{won}\t{total_wins[won]}")


if __name__ == "__main__":
    main()
