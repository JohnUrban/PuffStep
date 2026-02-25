#!/usr/bin/env python3
import sys, gzip, argparse, math

parser = argparse.ArgumentParser(
    description="""
Collapse consecutive bedGraph bins that share the same value in column 4.

Important:
- Only consecutive rows are considered
- Coordinates do NOT need to touch
- File must be sorted
"""
)

parser.add_argument("infile", help="Input bedGraph (or .gz)")
parser.add_argument("-o","--outfile", default="-",
                    help="Output file (default: stdout)")

parser.add_argument(
    "--tolerance",
    type=float,
    default=0.0,
    help="Treat values within this difference as equal (default: exact match)"
)

args = parser.parse_args()

# ---------- open ----------
def op(x):
    return gzip.open(x,"rt") if x.endswith(".gz") else open(x)

IN  = op(args.infile)
OUT = sys.stdout if args.outfile=="-" else open(args.outfile,"w")

# ---------- helper ----------
def same(a,b):
    if args.tolerance==0:
        return a==b
    return abs(a-b)<=args.tolerance

# ---------- collapse ----------
prev_chr=None
prev_start=None
prev_end=None
prev_val=None

for line in IN:
    if not line.strip():
        continue

    chrom,start,end,val = line.split()[:4]
    start=int(start)
    end=int(end)
    val=float(val)

    # first row
    if prev_chr is None:
        prev_chr,prev_start,prev_end,prev_val = chrom,start,end,val
        continue

    # same run?
    if chrom==prev_chr and same(val,prev_val):
        prev_end = end
        continue

    # flush previous
    OUT.write(f"{prev_chr}\t{prev_start}\t{prev_end}\t{prev_val}\n")

    # start new run
    prev_chr,prev_start,prev_end,prev_val = chrom,start,end,val

# final flush
if prev_chr is not None:
    OUT.write(f"{prev_chr}\t{prev_start}\t{prev_end}\t{prev_val}\n")

IN.close()
if OUT is not sys.stdout:
    OUT.close()

