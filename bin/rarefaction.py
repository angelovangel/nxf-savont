#!/usr/bin/env python3
"""
Rarefaction Curves — JSON output
-------------------------------------------
Usage:
    python rarefaction.py \
        --lineage  lineage-combined.tsv \
        --readstats d6306-A1_readstats.tsv d6306-A11_readstats.tsv \
        [--steps 50] [--iterations 20] [--max-depth 0] [--rarefaction rarefaction.json]

Requirements:
    pip install pandas numpy
"""

import argparse, json, os, sys
import numpy as np
import pandas as pd

METADATA_COLS = {"tax_id","species","genus","family","order","class","phylum","clade","superkingdom"}

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lineage",    required=True)
    p.add_argument("--readstats",  required=True, nargs="+")
    p.add_argument("--steps",      type=int, default=50)
    p.add_argument("--iterations", type=int, default=20)
    p.add_argument("--max-depth",  type=int, default=0)
    p.add_argument("--rarefaction",default="rarefaction.json")
    return p.parse_args()

def load_lineage(path):
    df = pd.read_csv(path, sep="\t")
    sample_cols = [c for c in df.columns if c not in METADATA_COLS]
    return df, sample_cols

def load_readstats(paths):
    total_reads = {}
    for p in paths:
        rs = pd.read_csv(p, sep="\t")
        for _, row in rs.iterrows():
            name = os.path.basename(str(row["file"])).split(".")[0]
            reads = int(row.get("num_seqs", row.get("reads", 0)))
            total_reads[name] = reads
    return total_reads

def build_count_table(df, sample_cols, total_reads):
    counts, missing = {}, []
    for s in sample_cols:
        if s in total_reads:
            counts[s] = np.round(df[s].values * total_reads[s]).astype(int)
        else:
            missing.append(s)
    if missing:
        print(f"[warn] No readstats for: {missing}", file=sys.stderr)
    if not counts:
        sys.exit("[error] No matching samples.")
    return counts

def rarefaction_numpy(counts, steps, n_iter, max_depth):
    result = {}
    for s, cv in counts.items():
        cap    = max_depth if max_depth > 0 else int(cv.sum())
        pool   = np.repeat(np.arange(len(cv)), cv)
        depths = np.unique(np.linspace(1, cap, steps).astype(int))
        rows   = []
        for d in depths:
            if d > len(pool):
                break
            richness = [len(np.unique(np.random.choice(pool, size=d, replace=False))) for _ in range(n_iter)]
            rows.append({"depth": int(d), "mean": float(np.mean(richness)), "std": float(np.std(richness))})
        result[s] = rows
    return result



def main():
    args = parse_args()
    np.random.seed(42)
    df, sample_cols = load_lineage(args.lineage)
    total_reads     = load_readstats(args.readstats)
    counts          = build_count_table(df, sample_cols, total_reads)
    print(f"[info] Samples: {list(counts.keys())}  |  Taxa: {len(df)}  |  Steps: {args.steps}  |  Iterations: {args.iterations}")
    series = rarefaction_numpy(counts, args.steps, args.iterations, args.max_depth)
    with open(args.rarefaction, "w") as f:
        json.dump(series, f)
    print(f"[done] JSON data -> {args.rarefaction}")

if __name__ == "__main__":
    main()
