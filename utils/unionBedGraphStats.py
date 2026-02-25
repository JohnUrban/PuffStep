#!/usr/bin/env python3
import sys
import argparse
import statistics
import gzip
import math

parser = argparse.ArgumentParser(
    description="Compute per-bin statistics across multiple bedGraph files."
)

parser.add_argument(
    "files",
    nargs="+",
    help="Input bedGraph files (can be .gz)"
)

parser.add_argument(
    "--check",
    action="store_true",
    help="Verify coordinates match across files"
)

parser.add_argument(
    "--skip-nan",
    action="store_true",
    help="Ignore NaN or '.' values instead of failing"
)

parser.add_argument(
    "--header",
    action="store_true",
    help="Return output with a header line"
)

args = parser.parse_args()

if args.header:
    print("chrom", "start", "end", "median", "mean", "min", "max", "mad", "iqr", "n", "stdev", "cv", sep="\t")

# ---------- file opener ----------
def openfile(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


handles = [openfile(f) for f in args.files]

line_num = 0

for lines in zip(*handles):
    line_num += 1

    rows = [l.rstrip().split() for l in lines]

    chrom, start, end = rows[0][0:3]

    # ---- coordinate check ----
    if args.check:
        for r in rows[1:]:
            if r[0] != chrom or r[1] != start or r[2] != end:
                sys.exit(
                    f"Coordinate mismatch at line {line_num}\n"
                    f"{rows[0]} vs {r}"
                )

    # ---- parse values ----
    vals = []
    for r in rows:
        v = r[3]
        if v in (".", "nan", "NaN"):
            if args.skip_nan:
                continue
            else:
                sys.exit(f"Missing value at line {line_num}")
        vals.append(float(v))

    n = len(vals)

    if n == 0:
        print(chrom, start, end, "NA","NA","NA","NA","NA","NA","NA", sep="\t")
        continue

    vals.sort()

    # ---- statistics ----
    median = vals[n//2] if n % 2 else (vals[n//2-1] + vals[n//2]) / 2
    mean = sum(vals)/n
    mn = vals[0]
    mx = vals[-1]

    # MAD
    mad = statistics.median(abs(v - median) for v in vals)

    # IQR
    q1 = statistics.median(vals[:n//2])
    q3 = statistics.median(vals[(n+1)//2:])
    iqr = q3 - q1

    # stdev + CV
    if n > 1:
        stdev = statistics.stdev(vals)
        cv = stdev / mean if mean != 0 else float("nan")
    else:
        stdev = 0.0
        cv = float("nan")

    print(
        chrom, start, end,
        median, mean, mn, mx,
        mad, iqr, n, stdev, cv,
        sep="\t"
    )


for h in handles:
    h.close()

