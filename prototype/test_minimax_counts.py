#!/usr/bin/env python3
"""Check minimax_counts against independent computations on a small graph.

Usage: python3 -I test_minimax_counts.py GRAPH.gfa PREFIX K [KMC_DUMP.txt]

GRAPH.gfa: the graph as `vg convert -f GRAPH.gbz` writes it (S and W lines).
PREFIX: the minimax_counts output prefix (needs PREFIX.dump.tsv from --dump-text and PREFIX.mmc).
KMC_DUMP.txt: `kmc_dump` of a KMC database counted from the same reads as the KFF (kmer, count).

genome_count is recomputed by walking every haplotype: the distinct graph positions of each canonical k-mer,
with gbwtgraph's rules (canonical = smaller 2-bit key; the position is the first base of the k-mer read in its
canonical orientation, as (node, strand, offset)). reads_count is compared with the KMC dump after
canonicalising its k-mers the same way. Prints numbers only.
"""
import collections
import struct
import sys

CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
COUNT_MAX = 65535


def encode(kmer):
    key = 0
    for c in kmer:
        key = (key << 2) | CODE[c]
    return key


def reverse_complement(key, k):
    rc = 0
    for _ in range(k):
        rc = (rc << 2) | (3 - (key & 3))
        key >>= 2
    return rc


def genome_counts(gfa, k):
    nodes, walks = {}, []
    with open(gfa) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "S":
                nodes[fields[1]] = fields[2].upper()
            elif fields[0] == "W":
                walk = fields[6]
                steps, i = [], 0
                while i < len(walk):
                    j = i + 1
                    while j < len(walk) and walk[j] not in "<>":
                        j += 1
                    steps.append((walk[i + 1:j], walk[i] == "<"))
                    i = j
                walks.append(steps)
    complement = str.maketrans("ACGT", "TGCA")
    positions = collections.defaultdict(set)
    for steps in walks:
        seq, pos = [], []  # per base: (node, is_reverse, offset) of that base on the walk's strand
        for node, rev in steps:
            s = nodes[node]
            if rev:
                s = s.translate(complement)[::-1]
            seq.append(s)
            pos += [(node, rev, o) for o in range(len(s))]
        seq = "".join(seq)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if any(c not in CODE for c in kmer):
                continue
            fwd = encode(kmer)
            rc = reverse_complement(fwd, k)
            if fwd < rc:
                key, p = fwd, pos[i]
            else:
                # Canonical is the reverse complement: it starts at the last base, read on the other strand.
                node, rev, off = pos[i + k - 1]
                key, p = rc, (node, not rev, len(nodes[node]) - 1 - off)
            positions[key].add(p)
    return {key: min(len(p), COUNT_MAX) for key, p in positions.items()}


def main():
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    gfa, prefix, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
    dump = {}
    with open(f"{prefix}.dump.tsv") as f:
        next(f)
        for line in f:
            key, g, r = map(int, line.split("\t"))
            dump[key] = (g, r)

    expected = genome_counts(gfa, k)
    print(f"graph_kmers\tprogram={len(dump)}\tbrute_force={len(expected)}")
    same_keys = dump.keys() == expected.keys()
    genome_mismatch = sum(1 for key in dump if key in expected and dump[key][0] != expected[key])
    print(f"genome_count\tsame_key_set={same_keys}\tmismatches={genome_mismatch}"
          f"\tgenome_count_1={sum(1 for g in expected.values() if g == 1)}\tmax={max(expected.values())}")
    failed = not same_keys or genome_mismatch > 0

    if len(sys.argv) == 5:
        kmc = {}
        with open(sys.argv[4]) as f:
            for line in f:
                kmer, count = line.split()
                key = encode(kmer)
                key = min(key, reverse_complement(key, k))
                kmc[key] = kmc.get(key, 0) + int(count)
        reads_mismatch = sum(1 for key, (g, r) in dump.items() if r != min(kmc.get(key, 0), COUNT_MAX))
        in_graph = sum(1 for key in kmc if key in dump)
        print(f"reads_count\tkmc_kmers={len(kmc)}\tkmc_kmers_in_graph={in_graph}"
              f"\tgraph_kmers_with_reads={sum(1 for g, r in dump.values() if r > 0)}\tmismatches={reads_mismatch}")
        failed |= reads_mismatch > 0

    # The binary counts file must hold the same table.
    with open(f"{prefix}.mmc", "rb") as f:
        magic, k_file, n = f.read(8), *struct.unpack("<QQ", f.read(16))
        keys = struct.unpack(f"<{n}Q", f.read(8 * n))
        packed = struct.unpack(f"<{n}I", f.read(4 * n))
        rest = f.read()
    file_ok = (magic == b"MMXCNT01" and k_file == k and n == len(dump) and not rest and list(keys) == sorted(dump)
               and all(dump[key] == (p >> 16, p & 0xFFFF) for key, p in zip(keys, packed)))
    print(f"counts_file\tmagic={magic.decode()}\tk={k_file}\tentries={n}\tbytes={24 + 12 * n}\tmatches_dump={file_ok}")
    failed |= not file_ok
    print("FAIL" if failed else "PASS")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
