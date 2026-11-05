#!/usr/bin/env python3
"""Tests for minimaxmer.py on synthetic sequences only.

Usage: python3 -I test_minimaxmer.py [--oracle PATH]
--oracle is a vg_oracle binary built from vg_oracle.cpp; without it the tests that
compare against gbwtgraph itself are skipped.
"""
import argparse
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import minimaxmer as mm  # noqa: E402

RNG = random.Random(20261020)
COMPLEMENT = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def random_seq(n, alphabet="ACGT"):
    return "".join(RNG.choice(alphabet) for _ in range(n))


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def oracle_lines(oracle, lines):
    out = subprocess.run([oracle], input="".join(lines), capture_output=True, text=True, check=True).stdout
    return out.splitlines()


def test_encoding_and_hash_against_gbwtgraph(oracle):
    kmers = ["A" * 31, "C" * 31, "ACGT" * 7 + "ACG", "TTGCA" * 6 + "G", random_seq(31), random_seq(31),
             "ACGTA", "acgta", random_seq(15)]
    rows = oracle_lines(oracle, [f"key {len(k)} {k}\n" for k in kmers])
    for kmer, row in zip(kmers, rows):
        fk, fh, rk, rh = map(int, row.split("\t"))
        k = len(kmer)
        assert mm.encode(kmer) == fk, kmer
        assert mm.wang_hash_64(fk) == fh, kmer
        assert mm.reverse_complement_code(fk, k) == rk, kmer
        assert mm.wang_hash_64(rk) == rh, kmer
    return f"{len(kmers)} k-mers: code, hash, reverse-complement code and hash equal gbwtgraph"


def test_rolling_encoding():
    seq = random_seq(500)
    for k in (5, 21, 31):
        for i, kmer in enumerate(mm.canonical_kmers(seq, k)):
            f = mm.encode(seq[i:i + k])
            r = mm.reverse_complement_code(f, k)
            want = (r, mm.wang_hash_64(r), True) if mm.wang_hash_64(r) < mm.wang_hash_64(f) else (f, mm.wang_hash_64(f), False)
            assert kmer == want, (k, i)
    return "rolling forward/reverse codes equal direct encoding (k = 5, 21, 31)"


def test_constant_counts_equal_gbwtgraph(oracle):
    cases = []
    for n in (60, 200, 1000, 5000):
        cases.append(random_seq(n))
    cases.append(random_seq(300) + "N" + random_seq(300) + "NNNN" + random_seq(300))   # invalid characters
    cases.append(random_seq(200).lower() + random_seq(200))                             # lowercase
    cases.append(("ACGTTGCA" * 100) + random_seq(200))                                   # tandem repeat: many ties
    cases.append(("A" * 300) + random_seq(100))                                          # homopolymer
    checked = 0
    for k, w in ((31, 20), (31, 10), (31, 50), (15, 5)):
        lines = [f"min {k} {w} {s}\n" for s in cases]
        out = oracle_lines(oracle, lines)
        results, current = [], []
        for row in out:
            if row == "end":
                results.append(current)
                current = []
            else:
                off, strand, key, h = row.split("\t")
                start = int(off) - (k - 1) if strand == "-" else int(off)
                current.append((start, strand == "-", int(key)))
        for seq, want in zip(cases, results):
            got, _ = mm.select_minimaxmers(seq, k, w)
            assert sorted(got) == sorted(want), (k, w, len(seq), len(got), len(want))
            checked += 1
    return f"{checked} sequence/(k, w) cases: selected offsets, strands and codes identical to gbwtgraph"


def test_density():
    seq = random_seq(1_000_000)
    lines = []
    for w in (10, 20, 50):
        selected, _ = mm.select_minimaxmers(seq, 31, w)
        density = len(selected) / len(seq)
        theory = 2 / (w + 1)
        assert abs(density / theory - 1) < 0.02, (w, density, theory)
        lines.append(f"w={w}: {density:.5f} vs 2/(w+1) = {theory:.5f} ({100 * (density / theory - 1):+.2f}%)")
    return "density on 1 Mbp random sequence, k=31: " + "; ".join(lines)


def test_genome_count_avoids_repetitive_kmer():
    seq = random_seq(3000)
    k, w = 31, 20
    plain, _ = mm.select_minimaxmers(seq, k, w)
    repetitive = {mm.count_key(code, k) for _, _, code in plain[::3]}  # pretend every third minimizer is repetitive
    counts = lambda key: 1000 if key in repetitive else 1
    got, _ = mm.select_minimaxmers(seq, k, w, genome_count=counts)
    chosen = {mm.count_key(code, k) for _, _, code in got}
    assert not (chosen & repetitive), "a high-genome-count k-mer was selected while a unique one was available"
    return f"{len(repetitive)} k-mers given genome count 1000 are never selected (k=31, w=20)"


def test_reads_count_breaks_ties():
    seq = random_seq(3000)
    k, w = 31, 20
    plain, _ = mm.select_minimaxmers(seq, k, w)
    plain_keys = {mm.count_key(code, k) for _, _, code in plain}
    kmers = mm.canonical_kmers(seq, k)
    # Non-minimizers at least 2w apart, so no window contains two of them.
    favoured, last = set(), -10 * w
    for off, kmer in enumerate(kmers):
        if kmer and mm.count_key(kmer[0], k) not in plain_keys and off - last >= 2 * w and len(favoured) < 20:
            favoured.add(mm.count_key(kmer[0], k))
            last = off
    got, _ = mm.select_minimaxmers(seq, k, w, reads_count=lambda key: 50 if key in favoured else 1)
    chosen = {mm.count_key(code, k) for _, _, code in got}
    assert favoured <= chosen, "a k-mer with a higher reads count lost against an equal genome count"
    return "20 k-mers that plain minimizers skip are all selected once their reads count is higher"


def test_reverse_complement_symmetry():
    for n in (500, 5000):
        seq = random_seq(n)
        a, _ = mm.select_minimaxmers(seq, 31, 20)
        b, _ = mm.select_minimaxmers(revcomp(seq), 31, 20)
        assert {code for _, _, code in a} == {code for _, _, code in b}, n
        assert sorted(off for off, _, _ in a) == sorted(n - 31 - off for off, _, _ in b), n
    return "a sequence and its reverse complement select the same canonical k-mers at mirrored offsets"


def test_edge_cases():
    k, w = 31, 20
    assert mm.select_minimaxmers(random_seq(k + w - 2), k, w) == ([], {})
    seq = random_seq(400)
    with_n = seq[:200] + "N" + seq[201:]
    got, _ = mm.select_minimaxmers(with_n, k, w)
    assert all(not (off <= 200 < off + k) for off, _, _ in got), "a k-mer containing N was selected"
    plain, _ = mm.select_minimaxmers(seq, k, w)
    absent = mm.count_key(plain[0][2], k)
    skip, _ = mm.select_minimaxmers(seq, k, w, genome_count=lambda c: 0 if c == absent else 1)
    keep, _ = mm.select_minimaxmers(seq, k, w, genome_count=lambda c: 0 if c == absent else 1, skip_absent=False)
    assert absent not in {mm.count_key(c, k) for _, _, c in skip} and absent in {mm.count_key(c, k) for _, _, c in keep}
    return "short sequence -> nothing; k-mers with N never selected; genome_count 0 skipped unless skip_absent=False"


def test_window_counts():
    k, w = 31, 20
    seq = random_seq(2000)
    selected, wins = mm.select_minimaxmers(seq, k, w)
    windows = len(seq) - (k + w - 1) + 1
    assert set(wins) == {off for off, _, _ in selected}
    assert sum(wins.values()) >= windows and min(wins.values()) >= 1
    return f"every selected k-mer wins >= 1 window; total wins {sum(wins.values())} for {windows} windows"


def repetitive_seq():
    """Random sequence with short tandem repeats (a k-mer recurs within one window, so windows have tied
    minima) and a dispersed repeat."""
    short, other, block = random_seq(7), random_seq(5), random_seq(300)
    return (random_seq(400) + short * 30 + random_seq(200) + other * 40 + random_seq(150) + block
            + random_seq(150) + block + random_seq(300))


def test_robust_winners_without_ties():
    k, w = 31, 20
    seq = random_seq(3000)
    plain, _ = mm.select_minimaxmers(seq, k, w)
    keys = [None if x is None else (1, mm.READS_COUNT_MAX, x[1] >> 32) for x in mm.canonical_kmers(seq, k)]
    nw = mm.robust_window_counts(seq, k, w, keys)
    windows = len(keys) - w + 1
    assert set(nw) == {off for off, _, _ in plain}, "without ties, robust winners must be the plain minimizers"
    assert sum(nw.values()) == windows, "every window must have exactly one winner"
    return f"no ties: winners = plain minimizers, {sum(nw.values())} wins for {windows} windows"


def test_robust_winners_with_ties():
    k, w = 31, 20
    seq = repetitive_seq()
    keys = [None if x is None else (1, mm.READS_COUNT_MAX, x[1] >> 32) for x in mm.canonical_kmers(seq, k)]
    nw = mm.robust_window_counts(seq, k, w, keys)
    windows = len(keys) - w + 1
    tied = sum(1 for s in range(windows) if keys[s:s + w].count(min(keys[s:s + w])) > 1)
    assert tied > 0, "the test sequence must have tied windows"
    extra = sum(nw.values()) - windows
    assert 0 <= extra, "a window lost its winner"
    rc = revcomp(seq)
    keys_rc = [None if x is None else (1, mm.READS_COUNT_MAX, x[1] >> 32) for x in mm.canonical_kmers(rc, k)]
    nw_rc = mm.robust_window_counts(rc, k, w, keys_rc)
    mirrored = {len(seq) - k - p: n for p, n in nw.items()}
    assert mirrored == dict(nw_rc), "window winners must mirror between a sequence and its reverse complement"
    return f"{tied} windows with tied minima; one winner per window (+{extra} for self-RC windows); winners mirror under RC"


def test_sparse_thinning():
    k, w = 31, 20
    for seq in (random_seq(4000), repetitive_seq()):
        for rspace in (0, 20, 100):
            accepted, nw = mm.select_sparse(seq, k, w, rspace)
            assert set(accepted) <= set(nw), "only window winners can be accepted"
            if rspace == 0:
                assert set(accepted) == set(nw), "rspace 0 must keep every window winner"
            else:
                assert all(b - a >= rspace for a, b in zip(accepted, accepted[1:])), "accepted seeds closer than rspace"
            rc_accepted, _ = mm.select_sparse(revcomp(seq), k, w, rspace)
            assert sorted(len(seq) - k - p for p in accepted) == rc_accepted, "thinning must mirror under RC"
    seq = random_seq(5000)
    counts = [len(mm.select_sparse(seq, k, w, r)[0]) for r in (0, 20, 100)]
    assert counts[0] > counts[1] > counts[2], counts
    return f"spacing respected, RC-mirrored, rspace 0 keeps all winners; 5 kbp seeds at rspace 0/20/100: {counts}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--oracle")
    args = parser.parse_args()
    tests = [test_rolling_encoding, test_density, test_genome_count_avoids_repetitive_kmer,
             test_reads_count_breaks_ties, test_reverse_complement_symmetry, test_edge_cases, test_window_counts,
             test_robust_winners_without_ties, test_robust_winners_with_ties, test_sparse_thinning]
    oracle_tests = [test_encoding_and_hash_against_gbwtgraph, test_constant_counts_equal_gbwtgraph]
    failed = 0
    for test in oracle_tests + tests:
        name = test.__name__
        if test in oracle_tests and not args.oracle:
            print(f"SKIP {name} (no --oracle)")
            continue
        try:
            detail = test(args.oracle) if test in oracle_tests else test()
            print(f"PASS {name}: {detail}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {name}: {e}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
