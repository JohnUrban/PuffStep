#!/usr/bin/env python3
import argparse, gzip, os, sys
from statistics import median
from collections import defaultdict

# --------------------------------------------------
# ARGUMENTS
# --------------------------------------------------
p = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    description="""
remove_zero_bins.py

Purpose
-------
Filter genomic bins from multiple bedGraph files based on zero-valued bins
observed in a set of denominator samples. This is typically used before
computing ratios (e.g. test/control coverage ratios) so that bins lacking
signal in controls are removed from ALL tracks consistently.

The program:

  1) Reads numerator and denominator bedGraph files simultaneously
  2) Determines which bins should be removed based ONLY on denominator values
  3) Removes those bins from every input file (numerator + denominator)
  4) Writes filtered copies of all files to --outdir
  5) Optionally outputs BED masks describing removed and/or kept bins
  6) Optionally prints chromosome-level summary statistics

Key assumptions
---------------
• All input files share identical bin coordinates and ordering
• bedGraphs are sorted and aligned
• Fourth column is numeric signal value

Typical use case
----------------
Remove bins where controls have zero signal so downstream ratios do not
contain infinities or unstable values.

Rule selection (REQUIRED — choose one)
---------------------------------------
--any        remove bin if ANY denominator sample == 0
--all        remove bin if ALL denominator samples == 0
--median     remove bin if median denominator value == 0
--mean       same behavior as --all (kept for semantic clarity)
--maxcount N remove bin if >N denominator samples are zero
--maxprop P  remove bin if proportion of zero denominators >P

Examples
--------
Remove bins if any control is zero:
  --any

Remove bins if majority of controls are zero:
  --maxprop 0.5

Remove bins if at least 3 controls are zero:
  --maxcount 2
"""
)

# -------------------------
# required inputs
# -------------------------
p.add_argument(
    "--num",
    nargs="+",
    required=True,
    help="""
Numerator bedGraph files.
These will be filtered using zero-bin rules determined from denominator files.
Multiple files allowed.
"""
)

p.add_argument(
    "--den",
    nargs="+",
    required=True,
    help="""
Denominator bedGraph files.
These determine which bins are removed.
Zero values in these files trigger removal depending on chosen rule.
"""
)

p.add_argument(
    "--outdir",
    required=True,
    help="""
Directory where filtered bedGraphs will be written.

Each output file keeps its original filename but gains suffix:
    .zeroBinsRemoved.bedGraph
"""
)

# -------------------------
# filtering rule
# -------------------------
g = p.add_mutually_exclusive_group(required=True)

g.add_argument(
    "--any",
    action="store_true",
    help="Remove bin if ANY denominator sample has value 0."
)

g.add_argument(
    "--median",
    action="store_true",
    help="Remove bin if median denominator value is 0."
)

g.add_argument(
    "--all",
    action="store_true",
    help="Remove bin only if ALL denominator samples are 0."
)

g.add_argument(
    "--mean",
    action="store_true",
    help="Same behavior as --all (included for naming convenience)."
)

g.add_argument(
    "--maxcount",
    type=int,
    help="Remove bin if more than N denominator samples are zero."
)

g.add_argument(
    "--maxprop",
    type=float,
    help="Remove bin if proportion of zero denominators exceeds P (0–1)."
)

# -------------------------
# safety / validation
# -------------------------
p.add_argument(
    "--check",
    action="store_true",
    help="""
Verify that all files share identical coordinates at each line.
Recommended unless you are certain files are perfectly aligned.
"""
)

# -------------------------
# optional mask outputs
# -------------------------
p.add_argument(
    "--removed-bed",
    help="""
Write BED file listing removed bins.

Output columns:
    chrom start end zeroCount zeroProp
"""
)

p.add_argument(
    "--keep-bed",
    help="""
Write BED file listing bins that were retained after filtering.

Output columns:
    chrom start end zeroCount zeroProp
"""
)

# -------------------------
# diagnostics
# -------------------------
p.add_argument(
    "--chrom-summary",
    action="store_true",
    help="""
Print per-chromosome statistics to stderr:

    chrom  kept  removed  prop_removed

Useful for detecting problematic chromosomes or samples.
"""
)

args = p.parse_args()



# --------------------------------------------------
# HELPERS
# --------------------------------------------------
def op(path):
    return gzip.open(path,"rt") if path.endswith(".gz") else open(path)

def outname(path):
    base = os.path.basename(path)
    if base.endswith(".gz"):
        base = base[:-3]
    return os.path.join(args.outdir, base.replace(".bedGraph",".zeroBinsRemoved.bedGraph"))

os.makedirs(args.outdir, exist_ok=True)

# --------------------------------------------------
# OPEN FILES
# --------------------------------------------------
den_handles = [op(f) for f in args.den]
num_handles = [op(f) for f in args.num]
all_handles = den_handles + num_handles

out_handles = [open(outname(f),"w") for f in args.den + args.num]

removed_bed = open(args.removed_bed,"w") if args.removed_bed else None
kept_bed    = open(args.keep_bed,"w") if args.keep_bed else None

# --------------------------------------------------
# RULE LOGIC
# --------------------------------------------------
def remove_bin(vals):
    zeros = [v==0 for v in vals]
    zcount = sum(zeros)
    n = len(vals)

    if args.any:
        return zcount>0, zcount, n
    if args.all or args.mean:
        return zcount==n, zcount, n
    if args.median:
        return median(vals)==0, zcount, n
    if args.maxcount is not None:
        return zcount>args.maxcount, zcount, n
    if args.maxprop is not None:
        return (zcount/n)>args.maxprop, zcount, n

# --------------------------------------------------
# CHROM STATS
# --------------------------------------------------
if args.chrom_summary:
    stats = defaultdict(lambda: {"kept":0,"removed":0})

# --------------------------------------------------
# MAIN LOOP
# --------------------------------------------------
line = 0
for rows in zip(*all_handles):
    line += 1
    rows = [r.rstrip().split() for r in rows]

    chrom,start,end = rows[0][:3]

    if args.check:
        for r in rows[1:]:
            if r[0]!=chrom or r[1]!=start or r[2]!=end:
                sys.exit(f"Coordinate mismatch line {line}")

    den_vals = [float(r[3]) for r in rows[:len(den_handles)]]

    remove, zcount, n = remove_bin(den_vals)
    prop = zcount/n if n else 0

    # ---------- BED OUTPUT ----------
    if remove and removed_bed:
        removed_bed.write(f"{chrom}\t{start}\t{end}\t{zcount}\t{prop:.5f}\n")

    if (not remove) and kept_bed:
        kept_bed.write(f"{chrom}\t{start}\t{end}\t{zcount}\t{prop:.5f}\n")

    # ---------- STATS ----------
    if args.chrom_summary:
        stats[chrom]["removed" if remove else "kept"] += 1

    # ---------- FILTER ----------
    if remove:
        continue

    for r,h in zip(rows,out_handles):
        h.write("\t".join(r)+"\n")

# --------------------------------------------------
# CLEANUP
# --------------------------------------------------
for h in all_handles + out_handles:
    h.close()

if removed_bed:
    removed_bed.close()
if kept_bed:
    kept_bed.close()

# --------------------------------------------------
# SUMMARY
# --------------------------------------------------
if args.chrom_summary:
    print("\nChromosome summary:", file=sys.stderr)
    print("chrom\tkept\tremoved\tprop_removed", file=sys.stderr)
    for c in sorted(stats):
        k = stats[c]["kept"]
        r = stats[c]["removed"]
        tot = k+r
        prop = r/tot if tot else 0
        print(f"{c}\t{k}\t{r}\t{prop:.4f}", file=sys.stderr)

