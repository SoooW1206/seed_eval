#!/usr/bin/env python3
"""Compare two `vg view -aj` outputs read by read, ignoring timing fields.

Usage: python3 -I compare_gam_json.py A.json B.json
A field is ignored if any component of its path is "time" or "time_used" (giraffe's
provenance and aligner timings). Reads are matched by name. Prints counts only.
"""
import json
import sys


def flat(value, path=()):
    if isinstance(value, dict):
        out = {}
        for k, v in value.items():
            out.update(flat(v, path + (k,)))
        return out
    return {path: value}


def load(path):
    reads = {}
    for line in open(path):
        r = json.loads(line)
        fields = {p: v for p, v in flat(r).items() if not ({"time", "time_used"} & set(p))}
        reads[r["name"]] = fields
    return reads


a, b = load(sys.argv[1]), load(sys.argv[2])
print(f"reads\t{len(a)}\t{len(b)}\nsame_names\t{a.keys() == b.keys()}")
diff_reads, diff_fields = 0, {}
for name in a.keys() & b.keys():
    keys = {p for p in a[name].keys() | b[name].keys() if a[name].get(p) != b[name].get(p)}
    diff_reads += bool(keys)
    for p in keys:
        diff_fields[".".join(map(str, p))] = diff_fields.get(".".join(map(str, p)), 0) + 1
compared = sum(len(f) for f in a.values())
print(f"fields_compared\t{compared}\nreads_with_differences\t{diff_reads}")
for field, n in sorted(diff_fields.items(), key=lambda t: -t[1])[:15]:
    print(f"differs\t{n}\t{field}")
print("IDENTICAL" if diff_reads == 0 and a.keys() == b.keys() else "DIFFERENT")
