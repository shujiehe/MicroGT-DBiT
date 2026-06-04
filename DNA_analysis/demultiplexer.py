#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Assign DNA reads to spatial coordinates by extracting B/A barcodes from R2
and mapping B+A to a whitelist.

Outputs:
  combine.fq
  demultiplexer.log
"""

import argparse
import os
import gzip
import sys
from collections import Counter

try:
    import regex as _rx
except Exception:
    sys.stderr.write(
        "[FATAL] The regex package is required. Install it with conda or pip.\n"
    )
    raise

DNA_ALPHABET = "ACGT"


def opengz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def iter_fastq(fh):
    while True:
        l1 = fh.readline()
        if not l1:
            break
        s = fh.readline()
        p = fh.readline()
        q = fh.readline()
        if not q:
            break
        yield l1.rstrip(), s.rstrip(), p.rstrip(), q.rstrip()


def parse_read_id(header_line):
    s = header_line.strip()
    if not s or s[0] != "@":
        return None
    return s[1:].split()[0]


def load_whitelist(path):
    wl = {}
    order = {}

    with open(path, "r") as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                sys.exit(
                    f"[ERROR] Invalid whitelist format at line {line_no}: {line}\n"
                    f"Expected: <B+A>\\t<x>\\t<y>"
                )

            bc = parts[0].strip().upper()
            x = parts[1].strip()
            y = parts[2].strip()

            if bc in wl:
                sys.stderr.write(
                    f"[WARN] Duplicate whitelist key at line {line_no}; keeping first: {bc}\n"
                )
                continue

            order[bc] = len(order)
            wl[bc] = f"{x}_{y}"

    return wl, order


_PATTERN_CACHE = {}


def get_fuzzy_pattern(L2, max_edits=1):
    key = (L2.upper(), int(max_edits))
    pat = _PATTERN_CACHE.get(key)

    if pat is None:
        pat = _rx.compile(r"(?:%s){e<=%d}" % (_rx.escape(L2.upper()), max_edits))
        _PATTERN_CACHE[key] = pat

    return pat


def _collect_L2_candidates(seq, L2, max_edits=1, compiled_pat=None, expect_start=40, expect_end=70):
    s = seq.upper()
    q = L2.upper()
    m = len(q)

    delta = 1 if max_edits <= 1 else max_edits
    min_span = m - delta
    max_span = m + delta

    def _dist(start, span):
        end = start + span
        return abs(start - expect_start) + abs(end - expect_end)

    cand_map = {}

    pos = s.find(q)
    while pos != -1:
        key = (pos, m)
        prev = cand_map.get(key)
        if prev is None or 0 < prev:
            cand_map[key] = 0
        pos = s.find(q, pos + 1)

    if max_edits > 0:
        pattern = compiled_pat or get_fuzzy_pattern(q, max_edits)

        for mo in pattern.finditer(s, overlapped=True):
            subs, ins, dels = mo.fuzzy_counts
            edits = subs + ins + dels

            if edits > max_edits:
                continue

            start = mo.start()
            span = mo.end() - mo.start()

            if span < min_span or span > max_span:
                continue

            key = (start, span)
            prev = cand_map.get(key)
            if prev is None or edits < prev:
                cand_map[key] = edits

    out = []
    for (start, span), edits in cand_map.items():
        out.append((edits, _dist(start, span), start, span))

    out.sort()
    return out


def rescue_1mm_ref_lazy(bc, wl_set, order_map):
    L = len(bc)
    candidates = []

    for i in range(L):
        orig = bc[i]
        for ch in DNA_ALPHABET:
            if ch == orig:
                continue

            v = bc[:i] + ch + bc[i + 1:]
            if v in wl_set:
                candidates.append(v)

    if not candidates:
        return None

    candidates.sort(key=lambda x: order_map.get(x, 10**12))
    return candidates[0]


def rescue_1mm_lazy(bc, wl, wl_set, order_map):
    ref = rescue_1mm_ref_lazy(bc, wl_set, order_map)
    if ref is None:
        return None
    return wl[ref]


def map_to_whitelist(bc, wl, wl_set, order_map, allow_mismatch=True):
    if bc in wl_set:
        return wl[bc], False

    if allow_mismatch:
        spatial = rescue_1mm_lazy(bc, wl, wl_set, order_map)
        if spatial is not None:
            return spatial, True

    return None, False


def extract_by_L2_with_WL(
    seqR2_raw,
    L2,
    wl_set,
    order_map,
    b_len=8,
    a_len=8,
    linker_mismatch=2,
    compiled_pat=None,
    allow_bc_rescue=True,
    expect_start=40,
    expect_end=70
):
    if not seqR2_raw:
        return None

    s = seqR2_raw.upper()

    cands = _collect_L2_candidates(
        s,
        L2,
        max_edits=linker_mismatch,
        compiled_pat=compiled_pat,
        expect_start=expect_start,
        expect_end=expect_end
    )

    if not cands:
        return None

    best_fallback = None
    best_exact = None
    rescue_todo = []

    for edits, dist, start, span in cands:
        if start < b_len:
            continue

        if (start + span + a_len) > len(s):
            continue

        B = s[start - b_len:start]
        A = s[start + span:start + span + a_len]
        bc = B + A
        score_base = (edits, dist, start, span)

        if best_fallback is None or score_base < best_fallback[0]:
            best_fallback = (score_base, B, A)

        if bc in wl_set:
            if best_exact is None or score_base < best_exact[0]:
                best_exact = (score_base, B, A)
        else:
            if allow_bc_rescue:
                rescue_todo.append((score_base, bc, B, A))

    if best_fallback is None and best_exact is None:
        return None

    if best_exact is not None:
        _, B, A = best_exact
        return B, A, "exact_wl"

    if allow_bc_rescue and rescue_todo:
        best_rescue = None

        for score_base, bc, B, A in rescue_todo:
            ref = rescue_1mm_ref_lazy(bc, wl_set, order_map)
            if ref is None:
                continue

            rescued_order = order_map.get(ref, 10**12)
            score = score_base + (rescued_order,)

            if best_rescue is None or score < best_rescue[0]:
                best_rescue = (score, B, A)

        if best_rescue is not None:
            _, B, A = best_rescue
            return B, A, "rescue_wl"

    _, B, A = best_fallback
    return B, A, "fallback"


def _read_one_fastq_block(fh):
    h = fh.readline()
    if not h:
        return False

    s = fh.readline()
    p = fh.readline()
    q = fh.readline()

    return bool(q)


def _count_remaining_blocks(fh, max_probe=1000):
    n = 0

    while n < max_probe and _read_one_fastq_block(fh):
        n += 1

    return n


def main():
    ap = argparse.ArgumentParser(
        description="Demultiplex DNA reads using L2 anchoring and whitelist mapping."
    )

    ap.add_argument("--R1", required=True)
    ap.add_argument("--R2", required=True)
    ap.add_argument("--whitelist", required=True, help="Whitelist: <B+A>\\t<x>\\t<y>")
    ap.add_argument("--outdir", required=True)

    ap.add_argument("--linker2", required=True)
    ap.add_argument("--a_len", type=int, default=8)
    ap.add_argument("--b_len", type=int, default=8)
    ap.add_argument("--min_r1_len", type=int, default=30)

    ap.add_argument(
        "--linker_mismatch",
        type=int,
        default=2,
        choices=[0, 1, 2],
        help="Maximum edit distance for L2 matching."
    )

    ap.add_argument(
        "--bc_max_mismatch",
        type=int,
        default=1,
        choices=[0, 1],
        help="Maximum Hamming distance for barcode rescue."
    )

    args = ap.parse_args()

    L2_real = args.linker2.upper()

    if len(L2_real) != 30:
        sys.stderr.write(
            f"[WARN] L2 length is not 30 bp; observed length = {len(L2_real)}.\n"
        )

    os.makedirs(args.outdir, exist_ok=True)

    out_fq = os.path.join(args.outdir, "combine.fq")
    log_fp = os.path.join(args.outdir, "demultiplexer.log")

    wl, order_map = load_whitelist(args.whitelist)
    wl_set = set(wl.keys())

    compiled_pat = (
        get_fuzzy_pattern(L2_real, args.linker_mismatch)
        if args.linker_mismatch > 0
        else None
    )

    total = 0
    ok_ext = 0
    ok_wl = 0
    d_short = 0
    d_linker = 0
    corrected_bc = 0
    not_in_wl_counts = Counter()
    pair_mismatch = 0

    pick_exact = 0
    pick_rescue = 0
    pick_fallback = 0

    with opengz(args.R1) as fr1, opengz(args.R2) as fr2, open(
        out_fq, "w", encoding="ascii", newline="\n"
    ) as fout:
        idx = 0

        for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(iter_fastq(fr1), iter_fastq(fr2)):
            id1 = parse_read_id(h1)
            id2 = parse_read_id(h2)

            if (id1 is None) or (id2 is None) or (id1 != id2):
                pair_mismatch += 1
                continue

            total += 1

            if len(s1) < args.min_r1_len:
                d_short += 1
                continue

            ex = extract_by_L2_with_WL(
                s2,
                L2_real,
                wl_set=wl_set,
                order_map=order_map,
                b_len=args.b_len,
                a_len=args.a_len,
                linker_mismatch=args.linker_mismatch,
                compiled_pat=compiled_pat,
                allow_bc_rescue=(args.bc_max_mismatch == 1),
                expect_start=40,
                expect_end=70
            )

            if ex is None:
                d_linker += 1
                continue

            B, A, tag = ex

            if tag == "exact_wl":
                pick_exact += 1
            elif tag == "rescue_wl":
                pick_rescue += 1
            else:
                pick_fallback += 1

            ok_ext += 1
            bc = B + A

            spatial, corrected = map_to_whitelist(
                bc,
                wl,
                wl_set,
                order_map,
                allow_mismatch=(args.bc_max_mismatch == 1)
            )

            if spatial is None:
                not_in_wl_counts[bc] += 1
                continue

            if corrected:
                corrected_bc += 1

            ok_wl += 1
            idx += 1

            name = f"@{idx}|:_:|{spatial}"
            fout.write(f"{name}\n{s1}\n+\n{q1}\n")

        extra_r1 = _count_remaining_blocks(fr1, max_probe=1000)
        extra_r2 = _count_remaining_blocks(fr2, max_probe=1000)
        pair_truncated = 1 if (extra_r1 > 0 or extra_r2 > 0) else 0

    top10_offwl = not_in_wl_counts.most_common(10)

    with open(log_fp, "w", encoding="utf-8") as lg:
        lg.write("[demultiplexer] DNA demultiplexing\n")
        lg.write(f"R1: {args.R1}\n")
        lg.write(f"R2: {args.R2}\n")
        lg.write(f"whitelist: {args.whitelist}\n")
        lg.write(f"L2: {L2_real}\n")
        lg.write(f"a_len={args.a_len}; b_len={args.b_len}\n")
        lg.write(f"min_r1_len={args.min_r1_len}\n")
        lg.write(f"linker_mismatch={args.linker_mismatch}\n")
        lg.write(f"bc_max_mismatch={args.bc_max_mismatch}\n\n")

        lg.write(f"Total pairs: {total}\n")
        lg.write(f"Output reads: {ok_wl}\n")
        lg.write(f"Pair-ID mismatched: {pair_mismatch}\n")
        lg.write(f"Dropped short R1: {d_short}\n")
        lg.write(f"Dropped no/invalid L2: {d_linker}\n")
        lg.write(f"Extracted by L2: {ok_ext}\n")
        lg.write(f"Extracted but not in whitelist: {ok_ext - ok_wl}\n")
        lg.write(f"Corrected barcodes: {corrected_bc}\n")
        lg.write(f"Pair truncated by zip: {pair_truncated} (R1_extra={extra_r1}, R2_extra={extra_r2})\n\n")

        lg.write(f"Pick exact WL: {pick_exact}\n")
        lg.write(f"Pick rescue WL: {pick_rescue}\n")
        lg.write(f"Pick fallback: {pick_fallback}\n\n")

        if top10_offwl:
            lg.write("Top10 off-whitelist barcodes:\n")
            for bc, ct in top10_offwl:
                lg.write(f"{bc}\t{ct}\n")


if __name__ == "__main__":
    main()