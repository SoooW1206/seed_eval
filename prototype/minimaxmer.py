#!/usr/bin/env python3
"""Prototype of minimaxmer selection (Yoshi Suzuki's himapper_rp), steps 1-2 only.

A minimaxmer is chosen like a minimizer, but each k-mer is ranked by the key
    (genome_count ascending, reads_count descending, hash ascending)
instead of by its hash alone: prefer k-mers that are rare in the genome, then
k-mers that are common in the reads. Step 1 selects the minimum of every window
of w consecutive k-mers; step 2 counts how many windows each selected k-mer won.
Yoshi's step 3 (uniformly sparse thinning) is not implemented.

Everything follows vg_minimax (gbwtgraph minimizer.h and minimax.h on the minimaxmer
branch), so the two can be compared read by read:
  * 2-bit encoding A=0 C=1 G=2 T=3 (either case), first base in the high bits;
  * hash = wang_hash_64 of the encoding; the reported orientation of a k-mer is
    the strand with the smaller hash (forward on ties);
  * counts are looked up with the smaller 2-bit code of the two strands, as in
    minimax_counts (count_key());
  * the key is (min(genome_count, 65534), 65535 - min(reads_count, 65535), hash >> 32),
    i.e. vg_minimax's packed order value; the high 32 bits of the hash order ties
    like the full hash;
  * a k-mer containing any other character is not a candidate;
  * every tied minimum of a window is reported, each read offset at most once.
With constant counts the selection is gbwtgraph's (w, k)-minimizers (up to ties in
the high 32 bits of the hash). As in Yoshi's code, k-mers with genome_count == 0 are
not candidates (skip_absent): like a k-mer with an N, such a k-mer never enters the
window buffer, the window still slides, and a window with no candidate selects nothing.

Command line (statistics only, never sequences):
    vg view -aj reads.gam | jq -r '[.name, .sequence] | @tsv' | python3 -I minimaxmer.py --k 31 --w 20
reads "name<TAB>sequence" lines, counts the k-mers of all reads for reads_count
(or uses a constant reads_count with --constant-counts), uses a constant
genome_count of 1, writes per-read selection counts (--per-read FILE: read name,
length, selected) and prints a summary.
"""
import argparse
import collections
import sys

MASK64 = (1 << 64) - 1
GENOME_COUNT_MAX = 0xFFFE  # as vg_minimax: one below the 16-bit maximum
READS_COUNT_MAX = 0xFFFF
READS_BIN = lambda r: 0 if r == 0 else 1 if r < 10 else 2  # --reads-mode bins: 0 / 1-9 / >=10
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


def count_key(code, k):
    """The code counts are stored under: the smaller 2-bit code of the two strands."""
    return min(code, reverse_complement_code(code, k))


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
    """Steps 1-2. genome_count(key) and reads_count(key) take count_key() of a k-mer
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
        c = count_key(code, k)
        g = genome_count(c) if genome_count else 1
        r = reads_count(c) if reads_count else 0
        keys.append(None if (skip_absent and g == 0) else
                    (min(g, GENOME_COUNT_MAX), READS_COUNT_MAX - min(r, READS_COUNT_MAX), h >> 32))

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


BALANCE = {"A": -1, "C": -1, "G": 1, "T": 1, "a": -1, "c": -1, "g": 1, "t": 1}


def window_orientation(seq, start, span):
    """Yoshi's D(W) = #G + #T - #A - #C over seq[start:start+span]: >0 means traverse right-to-left,
    <0 left-to-right. For D == 0 (even spans), compare the window with its reverse complement base by
    base (exact_window_orientation); 0 means the window is its own reverse complement."""
    d = sum(BALANCE.get(c, 0) for c in seq[start:start + span])
    if d:
        return 1 if d > 0 else -1
    for j in range(span):
        forward = PACK.get(seq[start + j], 0)
        reverse = 3 - PACK.get(seq[start + span - 1 - j], 0)
        if forward != reverse:
            return 1 if forward < reverse else -1
    return 0


def robust_window_counts(seq, k, w, keys):
    """Step 2 as in Yoshi's minimaxmer.c (count_canonical_robust_minima, default build): every window of w
    k-mers has one winner. A unique minimum wins. Among tied minima, D(W) orients the run of tied windows:
    D < 0 goes left to right, keeping the previous winner while it is still a minimum of the window and
    otherwise taking the rightmost minimum; D > 0 mirrors this right to left with the leftmost minimum; a
    tied window equal to its own reverse complement credits both extreme minima. keys[i] is None if k-mer i
    is not a candidate. Returns {position: windows won}."""
    plen, span = len(keys), k + w - 1
    nw = collections.Counter()
    if plen < w:
        return nw
    dq = collections.deque()
    positive = {"start": None, "leftmost": []}
    previous_class, forward_previous = 0, None   # classes: 0 reset, 1 unique, -1 negative, 2 positive

    def flush_positive(anchor):
        previous = anchor
        for j in reversed(range(len(positive["leftmost"]))):
            start, minimum = positive["start"] + j, positive["leftmost"][j]
            if previous is not None and start <= previous < start + w and keys[previous] == keys[minimum]:
                winner = previous
            else:
                winner = minimum
            nw[winner] += 1
            previous = winner
        positive["start"], positive["leftmost"] = None, []

    for i in range(plen):
        if keys[i] is not None:
            while dq and keys[i] < keys[dq[-1]]:       # equal keys stay: the front is the leftmost minimum
                dq.pop()
            dq.append(i)
        start = i - w + 1
        while dq and dq[0] < start:
            dq.popleft()
        if start < 0:
            continue
        if not dq:
            if positive["start"] is not None:
                flush_positive(None)
            previous_class, forward_previous = 0, None
            continue
        leftmost = rightmost = dq[0]
        for p in dq:
            if keys[p] != keys[leftmost]:
                break
            rightmost = p
        if rightmost == leftmost:                       # unique minimum
            if positive["start"] is not None:
                flush_positive(leftmost)
            nw[leftmost] += 1
            previous_class, forward_previous = 1, leftmost
            continue
        orientation = window_orientation(seq, start, span)
        if orientation < 0:
            if positive["start"] is not None:
                flush_positive(None)
            if previous_class not in (-1, 1):
                forward_previous = None
            if forward_previous is not None and start <= forward_previous <= i and keys[forward_previous] == keys[leftmost]:
                winner = forward_previous
            else:
                winner = rightmost
            nw[winner] += 1
            previous_class, forward_previous = -1, winner
        elif orientation > 0:
            if previous_class != 2:
                positive["start"], positive["leftmost"] = start, []
            positive["leftmost"].append(leftmost)
            previous_class, forward_previous = 2, None
        else:
            if positive["start"] is not None:
                flush_positive(None)
            nw[leftmost] += 1
            nw[rightmost] += 1
            previous_class, forward_previous = 0, None
    if positive["start"] is not None:
        flush_positive(None)
    return nw


def sparse_thinning(keys, nw, rspace):
    """Step 3 as in Yoshi's minimaxmer.c (canonical thinning, default build): candidates are the positions
    that won a window, sorted by (windows won descending, key ascending), i.e. (NW desc, genome_count asc,
    reads_count desc, hash asc). Within a group of equal (NW, key), positions are taken in pairs from both
    ends: both are accepted if at least rspace apart, otherwise neither, and the pair moves inward past
    blocked positions. An accepted position blocks every position closer than rspace (radius rspace - 1).
    Returns the accepted positions, sorted."""
    order = sorted(nw, key=lambda p: (-nw[p], keys[p], p))
    accepted, blocked = [], set()

    def accept(p):
        accepted.append(p)
        if rspace > 0:
            blocked.update(range(p - rspace + 1, p + rspace))

    begin = 0
    while begin < len(order):
        end = begin + 1
        while end < len(order) and nw[order[end]] == nw[order[begin]] and keys[order[end]] == keys[order[begin]]:
            end += 1
        group = order[begin:end]                       # already in position order within the group
        left, right = 0, len(group) - 1
        while left <= right:
            while left <= right and group[left] in blocked:
                left += 1
            if left > right:
                break
            while left < right and group[right] in blocked:
                right -= 1
            if left == right:
                accept(group[left])
                break
            if group[right] - group[left] >= rspace:
                accept(group[left])
                accept(group[right])
            left += 1
            right -= 1
        begin = end
    return sorted(accepted)


def select_sparse(seq, k, w, rspace, genome_count=None, reads_count=None):
    """Steps 1-3 with Yoshi's window winners and spacing. Same keys as select_minimaxmers (skip_absent).
    Returns (accepted positions, nw)."""
    keys = []
    for kmer in canonical_kmers(seq, k):
        if kmer is None:
            keys.append(None)
            continue
        code, h, _ = kmer
        c = count_key(code, k)
        g = genome_count(c) if genome_count else 1
        r = reads_count(c) if reads_count else 0
        keys.append(None if g == 0 else (min(g, GENOME_COUNT_MAX), READS_COUNT_MAX - min(r, READS_COUNT_MAX), h >> 32))
    nw = robust_window_counts(seq, k, w, keys)
    return sparse_thinning(keys, nw, rspace), nw


def load_mmc_counts(path, wanted, k, chunk=1 << 27):
    """{key: (genome_count, reads_count)} for the wanted keys found in a minimax_counts file.
    Reads the sorted keys sequentially in chunks instead of searching the whole file per key."""
    import numpy as np
    with open(path, "rb") as f:
        header = f.read(24)
    if header[:8] != b"MMXCNT01" or int.from_bytes(header[8:16], "little") != k:
        sys.exit(f"{path}: not a minimax_counts file with k={k}")
    n = int.from_bytes(header[16:24], "little")
    keys = np.memmap(path, dtype="<u8", mode="r", offset=24, shape=(n,))
    values = np.memmap(path, dtype="<u4", mode="r", offset=24 + 8 * n, shape=(n,))
    query = np.array(sorted(wanted), dtype=np.uint64)
    found = {}
    for start in range(0, n, chunk):
        block = np.asarray(keys[start:start + chunk])
        q = query[np.searchsorted(query, block[0], side="left"):np.searchsorted(query, block[-1], side="right")]
        idx = np.searchsorted(block, q)
        hit = idx < len(block)
        hit[hit] = block[idx[hit]] == q[hit]
        for key, i in zip(q[hit].tolist(), idx[hit].tolist()):
            v = int(values[start + i])
            found[key] = (v >> 16, v & 0xFFFF)
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--w", type=int, default=20)
    parser.add_argument("--constant-counts", action="store_true",
                        help="constant reads_count too: the selection is then exactly vg's minimizers")
    parser.add_argument("--counts-mmc", help="genome_count and reads_count from this minimax_counts file "
                                             "(k-mers not in it have genome_count 0)")
    parser.add_argument("--reads-mode", choices=["exact", "bins", "constant"], default="exact",
                        help="with --counts-mmc: use reads_count as is (exact), only its bin 0 / 1-9 / >=10 "
                             "(bins; ties within a bin go to the hash), or not at all (constant: genome_count only)")
    parser.add_argument("--rspace", type=int, help="step 3 as in Yoshi's code: window winners by D-oriented "
                        "robust winnowing, then uniformly sparse thinning with this spacing (k-mer positions)")
    parser.add_argument("--per-read", help="write read name, length and number of selected k-mers here")
    args = parser.parse_args()
    if args.constant_counts and args.counts_mmc:
        sys.exit("--constant-counts and --counts-mmc exclude each other")
    if args.reads_mode != "exact" and not args.counts_mmc:
        sys.exit("--reads-mode needs --counts-mmc")

    names, sequences = [], []
    for line in sys.stdin:
        if line.strip():
            name, seq = line.rstrip("\n").split("\t")
            names.append(name)
            sequences.append(seq)
    if len(set(names)) != len(names):
        sys.exit("duplicate read names in the input")
    reads = collections.Counter()
    genome_count = None
    if args.counts_mmc:
        wanted = {count_key(kmer[0], args.k) for seq in sequences for kmer in canonical_kmers(seq, args.k) if kmer}
        counts = load_mmc_counts(args.counts_mmc, wanted, args.k)
        genome_count = lambda key: counts.get(key, (0, 0))[0]
        reads_count = {"exact": lambda key: counts.get(key, (0, 0))[1],
                       "bins": lambda key: READS_BIN(counts.get(key, (0, 0))[1]),
                       "constant": None}[args.reads_mode]
        print(f"distinct_kmers_in_reads\t{len(wanted)}\nkmers_in_counts_file\t{len(counts)}")
    elif not args.constant_counts:
        for seq in sequences:
            reads.update(count_key(kmer[0], args.k) for kmer in canonical_kmers(seq, args.k) if kmer is not None)
        reads_count = reads.__getitem__
    else:
        reads_count = None
    per_read, total_bp, total_selected, total_wins = [], 0, 0, collections.Counter()
    for name, seq in zip(names, sequences):
        if args.rspace is None:
            selected, wins = select_minimaxmers(seq, args.k, args.w, genome_count=genome_count, reads_count=reads_count)
            won = wins.values()
        else:
            selected, nw = select_sparse(seq, args.k, args.w, args.rspace, genome_count=genome_count, reads_count=reads_count)
            won = [nw[p] for p in selected]
        per_read.append((name, len(seq), len(selected)))
        total_bp += len(seq)
        total_selected += len(selected)
        total_wins.update(won)
    if args.per_read:
        with open(args.per_read, "w") as f:
            f.write("read\tlength\tselected\n")
            f.writelines(f"{name}\t{n}\t{s}\n" for name, n, s in per_read)
    source = ("constant" if args.constant_counts else f"counts file {args.counts_mmc}, reads {args.reads_mode}"
              if args.counts_mmc else "reads_count from these reads, genome_count 1")
    print(f"counts\t{source}")
    if args.rspace is not None:
        print(f"step3\trobust window winners, sparse thinning rspace={args.rspace}")
    print(f"reads\t{len(sequences)}\ntotal_bp\t{total_bp}\nselected\t{total_selected}")
    print(f"density\t{total_selected / total_bp:.6g}\ntheory_2/(w+1)\t{2 / (args.w + 1):.6g}")
    if not args.constant_counts and not args.counts_mmc:
        print(f"distinct_kmers_in_reads\t{len(reads)}\nkmers_seen_once\t{sum(1 for c in reads.values() if c == 1)}")
    print("windows_won\tselected_kmers")
    for won in sorted(total_wins):
        print(f"{won}\t{total_wins[won]}")


if __name__ == "__main__":
    main()
