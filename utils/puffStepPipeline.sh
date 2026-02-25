#!/bin/bash
set -euo pipefail



MANIFEST=${1}


MAIN=$( dirname $( dirname ${0} ) )
UTILS=${MAIN}/utils
PUFFSTEP=${MAIN}/puffStep.py
UNIONBDGSTATS=${UTILS}/unionBedGraphStats.py
RMZEROBINS=${UTILS}/removeZeroBins.py
COLLAPSE=${UTILS}/collapseBedGraphRuns.py


STEP1DIR=01-mednorm
STEP2DIR=02-unionStats
STEP3DIR=03-removeZeroBins
STEP4DIR=04-chromMedRatioNorm-RCN
STEP5DIR=05-medianSmoothed-RCN
STEP6DIR=06-HMM
STEP7DIR=07-collapsedHMM
STEP8DIR=08-aboveBackground
STEP9DIR=09-summitStates
STEP10DIR=10-summitBins

## TO OPTIONALLY ADD:
STEP11DIR=11-medianSmoothed-log2RCN
STEP12DIR=12-HMM-log2RCN
STEP13DIR=13-collapsedHMM-log2RCN
STEP14DIR=14-aboveBackground-log2RCN
STEP15DIR=15-summitStates-log2RCN
STEP16DIR=16-summitBins-log2RCN


## Pseudocount to add when performing median ratio normalization. 
## Default is 0, but can be set to a small value (e.g. 1) if you want to avoid division by zero or log2(0) issues in later steps. 
## However, adding a pseudocount can also affect the results, especially for bins with low counts, so use with caution and consider the implications for your specific dataset and analysis goals. 
## Morever, since we do zero bin removal from the denominator before this step, a pseudocount of 0 should not cause any issues with division by zero in the RCN calculation.
PSEUDO=0 

## Smoothing halfwidth for median smoothing step. This determines how many bins on either side of the center bin are included in the median calculation. 
##  For example, a halfwidth of 3 means that the median will be calculated using the center bin and the 3 bins on either side (total of 7 bins). 
##  Adjust this value based on the bin size of your data and the desired level of smoothing. 
##  For larger bins, you may want a smaller halfwidth, while for smaller bins, a larger halfwidth may be appropriate to achieve the desired smoothing effect.
SMOOTHING_HALFWIDTH=3 ## For 5 kb bins, this corresponds to a ~30 kb window (3 bins on either side + the center bin). Adjust as needed for different bin sizes or desired smoothing levels.

# PARALLEL_OPTS="--keep-order --halt soon,fail=1 --jobs 0"
PARALLEL_OPTS=(
  --keep-order
  --halt soon,fail=1
  --jobs 0
)


####################################################################################################################
####################################################################################################################
## STEP 1 MEDIAN NORM ALL SAMPLES
## Read in manifest - group sample 
## As each line is read in, perform median normalization and rename with group header. 
## Also build up a sorted list of unique group numbers to be used in later steps.
####################################################################################################################
####################################################################################################################

######################################################################
## STEP 1 — MEDIAN NORMALIZE ALL SAMPLES (PARALLEL SAFE)
######################################################################

mkdir -p "$STEP1DIR"
shopt -s nullglob

declare -A GROUP_SEEN=()
MANIFEST_LINES=()

# --------------------------------------------------
# PASS 1 — parse + validate manifest
# --------------------------------------------------
while read -r GROUP INBDG; do

    [[ -z "$GROUP" ]] && continue

    if ! [[ "$GROUP" =~ ^[0-9]+$ ]]; then
        echo "ERROR: Group '$GROUP' is not an integer" >&2
        exit 1
    fi

    [[ -f "$INBDG" ]] || {
        echo "ERROR: File not found: $INBDG" >&2
        exit 1
    }

    GROUP_SEEN["$GROUP"]=1
    MANIFEST_LINES+=( "$GROUP"$'\t'"$INBDG" )

done < "$MANIFEST"


# --------------------------------------------------
# Build sorted group list
# --------------------------------------------------
mapfile -t GROUP_LIST < <(
    printf "%s\n" "${!GROUP_SEEN[@]}" | sort -n
)


# --------------------------------------------------
# Validate numbering
# --------------------------------------------------
if [[ "${GROUP_LIST[0]}" != "1" ]]; then
    echo "ERROR: Group numbering must start at 1 (found ${GROUP_LIST[0]})" >&2
    exit 1
fi

expected=1
for g in "${GROUP_LIST[@]}"; do
    if [[ "$g" -ne "$expected" ]]; then
        echo "ERROR: Missing group number $expected (found $g instead)" >&2
        exit 1
    fi
    ((expected++))
done


# --------------------------------------------------
# PRINT GROUP SUMMARY  ← INSERT HERE
# --------------------------------------------------
# echo "Detected groups: ${GROUP_LIST[*]}" >&2
# echo "Detected groups (${#GROUP_LIST[@]}): ${GROUP_LIST[*]}" >&2
echo "[INFO] Detected groups (${#GROUP_LIST[@]}): ${GROUP_LIST[*]}" >&2



# --------------------------------------------------
# PASS 2 — PARALLEL NORMALIZATION
# --------------------------------------------------
export PUFFSTEP STEP1DIR

printf "%s\n" "${MANIFEST_LINES[@]}" | \
parallel "${PARALLEL_OPTS[@]}" --colsep '\t' '
    GROUP=$(cut -f1 <<< "{}")
    INBDG=$(cut -f2 <<< "{}")

    BASE=$(basename "$INBDG" .bedGraph)

    echo "Normalizing: group $GROUP — $BASE" >&2

    "$PUFFSTEP" normalize --protocol1 -l "$INBDG" \
        | LC_ALL=C sortBed -i - \
        > "$STEP1DIR/${GROUP}.${BASE}.medNorm.bedGraph"
'

shopt -u nullglob


## SIMPLER UN-PARALLEL CODE:
# mkdir -p "$STEP1DIR"

# declare -A GROUP_SEEN=()

# while read -r GROUP INBDG; do

#     # skip empty lines
#     [[ -z "$GROUP" ]] && continue

#     # validate integer group
#     if ! [[ "$GROUP" =~ ^[0-9]+$ ]]; then
#         echo "ERROR: Group '$GROUP' is not an integer" >&2
#         exit 1
#     fi

#     GROUP_SEEN["$GROUP"]=1

#     INBASE=$(basename "$INBDG" .bedGraph)

#     "$PUFFSTEP" normalize --protocol1 -l "$INBDG" | sortBed -i - > "$STEP1DIR/${GROUP}.${INBASE}.medNorm.bedGraph"

# done < "$MANIFEST"


# # --------------------------------------------------
# # Build sorted unique group list
# # --------------------------------------------------

# mapfile -t GROUP_LIST < <(
#     printf "%s\n" "${!GROUP_SEEN[@]}" | sort -n
# )

# # --------------------------------------------------
# # Validate group numbering
# # --------------------------------------------------

# if [[ "${GROUP_LIST[0]}" != "1" ]]; then
#     echo "ERROR: Group numbering must start at 1 (found ${GROUP_LIST[0]})" >&2
#     exit 1
# fi

# # ensure consecutive numbering
# expected=1
# for g in "${GROUP_LIST[@]}"; do
#     if [[ "$g" -ne "$expected" ]]; then
#         echo "ERROR: Missing group number $expected (found $g instead)" >&2
#         exit 1
#     fi
#     ((expected++))
# done


# # --------------------------------------------------
# # Optional diagnostic print
# # --------------------------------------------------
# echo "Detected groups: ${GROUP_LIST[*]}" >&2




####################################################################################################################
####################################################################################################################
## STEP 2 — UNION BEDGRAPH STATS + PER-BIN MEDIAN
## Take all median normalized bedGraphs and union them together to get per bin stats and per bin median across samples.
####################################################################################################################
####################################################################################################################


######################################################################
## STEP 2 — UNION BEDGRAPH STATS + PER-BIN MEDIAN (PARALLEL)
######################################################################

mkdir -p "$STEP2DIR"
shopt -s nullglob

GROUP_JOBS=()

# --------------------------------------------------
# validation phase — resolve files per group
# --------------------------------------------------
for N in "${GROUP_LIST[@]}"; do

    files=( "$STEP1DIR"/"$N"*.medNorm.bedGraph )

    if (( ${#files[@]} == 0 )); then
        echo "ERROR: No files found for group $N in $STEP1DIR" >&2
        exit 1
    fi

    # store job spec: group + files
    GROUP_JOBS+=( "$N:::${files[*]}" )

done


# --------------------------------------------------
# export vars for parallel
# --------------------------------------------------
export UNIONSTATS STEP2DIR


# --------------------------------------------------
# parallel execution
# --------------------------------------------------
printf "%s\n" "${GROUP_JOBS[@]}" | \
parallel "${PARALLEL_OPTS[@]}" --colsep ':::' '

    N={1}
    FILES=({2})

    echo "Processing group $N (${#FILES[@]} files)" >&2

    outstats="$STEP2DIR/${N}.TPM.medNorm.unionBedGraphStats.txt"
    outbed="$STEP2DIR/${N}.TPM.medNorm.stage_${N}_median.bedGraph"

    "$UNIONSTATS" --check --skip-nan "${FILES[@]}" > "$outstats"

    awk "BEGIN{OFS=\"\t\"} {print \$1,\$2,\$3,\$4}" "$outstats" \
        | LC_ALL=C sortBed -i - \
        > "$outbed"
'

shopt -u nullglob





## simpler code: Unparallel version that processes groups sequentially but still uses unionBedGraphStats for each group:
# mkdir -p "$STEP2DIR"

# shopt -s nullglob   # unmatched globs expand to empty list

# for N in "${GROUP_LIST[@]}"; do

#     files=( "$STEP1DIR"/"$N"*".medNorm.bedGraph" )

#     # safety check
#     if (( ${#files[@]} == 0 )); then
#         echo "ERROR: No files found for group $N in $STEP1DIR" >&2
#         exit 1
#     fi

#     echo "Processing group $N (${#files[@]} files)" >&2

#     outstats="$STEP2DIR/${N}.TPM.medNorm.unionBedGraphStats.txt"
#     outbed="$STEP2DIR/${N}.TPM.medNorm.stage_${N}_median.bedGraph"

#     "$UNIONSTATS" --check --skip-nan "${files[@]}" > "$outstats"

#     # column 4 = median from unionstats output
#     awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' "$outstats" | sortBed -i - > "$outbed"

# done

# shopt -u nullglob




####################################################################################################################
####################################################################################################################
## STEP 3 — REMOVE ZERO BINS
## Take in the median bedGraph files and remove bins with zero values. 
## Output both removed and kept bins as separate bed files for use in later steps or QC. 
## Also output a chrom summary file with number of bins removed/kept per chrom and median value of removed bins per chrom.
####################################################################################################################
####################################################################################################################


mkdir -p "$STEP3DIR"

shopt -s nullglob

# denominator = files starting with "1."
den_files=( "$STEP2DIR"/1.*.bedGraph )

# numerator = all other bedGraphs
num_files=()
for f in "$STEP2DIR"/*.bedGraph; do
    base=$(basename "$f")
    [[ $base == 1.* ]] && continue
    num_files+=( "$f" )
done


# safety checks
if (( ${#den_files[@]} == 0 )); then
    echo "ERROR: no denominator files (1.*.bedGraph) found" >&2
    exit 1
fi

if (( ${#num_files[@]} == 0 )); then
    echo "ERROR: no numerator files found" >&2
    exit 1
fi


# run tool
"$ZEROBINRM" \
    --median \
    --outdir "$STEP3DIR" \
    --removed-bed "$STEP3DIR/removed.bed" \
    --keep-bed "$STEP3DIR/kept.bed" \
    --chrom-summary \
    --num "${num_files[@]}" \
    --den "${den_files[@]}" \
    > "$STEP3DIR/chrom-summary.txt" 2>&1




####################################################################################################################
####################################################################################################################
## STEP 4 - PERFORM CHROMOSOME-SPECIFIC MEDIAN RATIO NORMALIZATION USING GROUP 1 AS DENOMINATOR
## Take in the zero-bin-removed bedGraphs and perform chromosome-specific median ratio normalization using group 1 as the denominator. 
## Output normalized bedGraphs.
## Strict version: exactly ONE test file and ONE control file allowed
####################################################################################################################
####################################################################################################################


######################################################################
## STEP 4 — CHROMOSOME-SPECIFIC MEDIAN RATIO NORMALIZATION
## Parallel version — strict single test + single control
######################################################################

mkdir -p "$STEP4DIR"
shopt -s nullglob

# --------------------------------------------------
# locate CONTROL file (must be exactly one)
# --------------------------------------------------
CONTROL_FILES=( "$D"/1.*.bedGraph )

if (( ${#CONTROL_FILES[@]} == 0 )); then
    echo "ERROR: No control file found matching 1.*.bedGraph in $D" >&2
    exit 1
elif (( ${#CONTROL_FILES[@]} > 1 )); then
    echo "ERROR: Multiple control files found:" >&2
    printf '  %s\n' "${CONTROL_FILES[@]}" >&2
    exit 1
fi

CONTROL="${CONTROL_FILES[0]}"

echo "Control file:"
echo "  $CONTROL"
echo >&2

# --------------------------------------------------
# resolve all test files first (validation phase)
# --------------------------------------------------
TEST_LIST=()

for N in "${GROUP_LIST[@]:1}"; do

    TEST_FILES=( "$D"/"$N".*.bedGraph )

    if (( ${#TEST_FILES[@]} == 0 )); then
        echo "ERROR: No test file found for group $N" >&2
        exit 1
    elif (( ${#TEST_FILES[@]} > 1 )); then
        echo "ERROR: Multiple test files found for group $N:" >&2
        printf '  %s\n' "${TEST_FILES[@]}" >&2
        exit 1
    fi

    TEST_LIST+=( "${TEST_FILES[0]}" )

done

# --------------------------------------------------
# parallel execution
# --------------------------------------------------
export PUFFSTEP STEP4DIR CONTROL PSEUDO

printf "%s\n" "${TEST_LIST[@]}" | parallel "${PARALLEL_OPTS[@]}" \
    --jobs 0 \
    --halt soon,fail=1 \
    '

    TEST="{}"
    BASE=$(basename "$TEST" .bedGraph)

    echo "RCN normalization: $BASE" >&2

    "$PUFFSTEP" normalize \
        --protocol31 \
        -l "$TEST" \
        -e "$CONTROL" \
        --pseudo "$PSEUDO" \
    | LC_ALL=C sortBed -i - \
    > "$STEP4DIR/${BASE}.RCN_pseudo${PSEUDO}.bedGraph"
    '

shopt -u nullglob



## simpler / not parallel code:

# mkdir -p "$STEP4DIR"
# shopt -s nullglob

# # --------------------------------------------------
# # locate CONTROL file (must be exactly one)
# # pattern 1.* ensures we match "1." but NOT 10.*,100.*,etc.
# # --------------------------------------------------
# CONTROL_FILES=( "$D"/1.*.bedGraph )

# if (( ${#CONTROL_FILES[@]} == 0 )); then
#     echo "ERROR: No control file found matching 1.*.bedGraph in $D" >&2
#     exit 1
# elif (( ${#CONTROL_FILES[@]} > 1 )); then
#     echo "ERROR: Multiple control files found:" >&2
#     printf '  %s\n' "${CONTROL_FILES[@]}" >&2
#     exit 1
# fi

# CONTROL="${CONTROL_FILES[0]}"

# echo "Control file:"
# echo "  $CONTROL"
# echo >&2


# # --------------------------------------------------
# # process test groups (2..N)
# # --------------------------------------------------
# for N in "${GROUP_LIST[@]:1}"; do

#     TEST_FILES=( "$D"/"$N".*.bedGraph )

#     if (( ${#TEST_FILES[@]} == 0 )); then
#         echo "ERROR: No test file found for group $N" >&2
#         exit 1
#     elif (( ${#TEST_FILES[@]} > 1 )); then
#         echo "ERROR: Multiple test files found for group $N:" >&2
#         printf '  %s\n' "${TEST_FILES[@]}" >&2
#         exit 1
#     fi

#     TEST="${TEST_FILES[0]}"
#     BASE=$(basename "$TEST" .bedGraph)

#     echo "--------------------------------------------------"
#     echo "RCN normalization — group $N"
#     echo "Test     : $TEST"
#     echo "Control  : $CONTROL"
#     echo "Pseudo   : $PSEUDO"

#     "$PUFFSTEP" normalize \
#         --protocol31 \
#         -l "$TEST" \
#         -e "$CONTROL" \
#         --pseudo "$PSEUDO" | sortBed -i - \
#         > "$STEP4DIR/${BASE}.RCN_pseudo${PSEUDO}.bedGraph"

#     echo
# done

# shopt -u nullglob



####################################################################################################################
####################################################################################################################
## STEP 5 - PERFORM MEDIAN SMOOTHING
## Take in the chromosome-specific median ratio normalized bedGraphs and perform median smoothing. 
## Output smoothed bedGraphs.
####################################################################################################################
####################################################################################################################

######################################################################
## STEP 5 — MEDIAN SMOOTHING (PARALLEL + DETERMINISTIC SORT)
######################################################################

mkdir -p "$STEP5DIR"
shopt -s nullglob

files=( "$STEP4DIR"/*.bedGraph )

if (( ${#files[@]} == 0 )); then
    echo "ERROR: No input bedGraph files found in $STEP4DIR" >&2
    exit 1
fi

# --------------------------------------------------
# export variables so GNU parallel can see them
# --------------------------------------------------
export PUFFSTEP STEP5DIR SMOOTHING_HALFWIDTH

# --------------------------------------------------
# run smoothing in parallel
# --------------------------------------------------
printf "%s\n" "${files[@]}" | parallel "${PARALLEL_OPTS[@]}" '
    RCN="{}"
    BASE=$(basename "$RCN" .bedGraph)

    echo "Smoothing: $BASE" >&2

    "$PUFFSTEP" normalize \
        --protocol32 \
        --halfwidth "$SMOOTHING_HALFWIDTH" \
        -l "$RCN" \
    | LC_ALL=C sortBed -i - \
    > "$STEP5DIR/${BASE}.medSmoothedHW${SMOOTHING_HALFWIDTH}.bedGraph"
'

shopt -u nullglob

## Simpler code:
# mkdir -p "$STEP5DIR"

# for RCN in "$STEP4DIR"/*bedGraph ; do 
#   BASE=$( basename "${RCN}" .bedGraph ) ; 
#   ${PUFFSTEP} normalize --protocol32 --halfwidth "${SMOOTHING_HALFWIDTH}" -l "${RCN}" | sortBed -i - > "$STEP5DIR/${BASE}.medSmoothedHW${SMOOTHING_HALFWIDTH}.bedGraph"
# done


####################################################################################################################
####################################################################################################################
## STEP 6 - HMM PUFFCN STATE CALLING
## Take in the median smoothed bedGraphs and run HMM to get statepath bedGraphs
## Output statepath bedGraphs.
####################################################################################################################
####################################################################################################################

# ${PUFFSTEP} puffcn -i ${F} --leave_special_state ${LEAVE_SPECIAL} --leave_other ${LEAVE_OTHER} --exp_decay --mu 1,2,4,8,16,32,64,128 --sigma 0.25,0.5,1,2,4,6,8,24 > ${BASEF}.expdecay.states.bedGraph
# ${PUFFSTEP} puffcn -i ${F} --leave_special_state ${LEAVE_SPECIAL} --leave_other ${LEAVE_OTHER} --exp_decay --mu 1,2,4,8,16,32,64,128 --sigma 0.25,0.5,1,2,4,8,16,32 > ${BASEF}.expdecay.states.bedGraph ;

mkdir -p "$STEP6DIR"

# -------------------------
# PARAMETERS
# -------------------------
EMISSION_MEANS="1,2,4,8,16,32,64,128"
EMISSION_SIGMAS="0.25,0.5,1,2,4,6,8,24"

BIN_SIZE=5000
EXPECTED_SPECIAL_LENGTH=1000000
EXPECTED_OTHER_LENGTH=75000

# compute transition probabilities (deterministic, no awk subshell parsing issues)
LEAVE_SPECIAL=$(awk -v b="$BIN_SIZE" -v e="$EXPECTED_SPECIAL_LENGTH" 'BEGIN{print b/e}')
LEAVE_OTHER=$(awk -v b="$BIN_SIZE" -v e="$EXPECTED_OTHER_LENGTH" 'BEGIN{print b/e}')

echo "Parameters:"
echo "  LEAVE_SPECIAL = $LEAVE_SPECIAL"
echo "  LEAVE_OTHER   = $LEAVE_OTHER"
echo

shopt -s nullglob
files=( "$STEP5DIR"/*.bedGraph )

if (( ${#files[@]} == 0 )); then
    echo "ERROR: no Step5 bedGraphs found in $STEP5DIR" >&2
    exit 1
fi

export PUFFSTEP EMISSION_MEANS EMISSION_SIGMAS LEAVE_SPECIAL LEAVE_OTHER

# --------------------------------------------------
# PARALLEL EXECUTION
# --------------------------------------------------
printf "%s\n" "${files[@]}" | parallel "${PARALLEL_OPTS[@]}" '

    F="{}"
    BASE=$(basename "$F" .bedGraph)

    echo "Calling states: $BASE" >&2

    "$PUFFSTEP" puffcn \
        -i "$F" \
        --leave_special_state "$LEAVE_SPECIAL" \
        --leave_other "$LEAVE_OTHER" \
        --exp_decay \
        --mu "$EMISSION_MEANS" \
        --sigma "$EMISSION_SIGMAS" \
    | LC_ALL=C sortBed -i - \
    > "$STEP6DIR/${BASE}.expdecay.states.bedGraph"
'

shopt -u nullglob





## Simpler code:

# EMISSION_MEANS="1,2,4,8,16,32,64,128"
# EMISSION_SIGMAS="0.25,0.5,1,2,4,6,8,24"
# BIN_SIZE=5000 ; 
# EXPECTED_SPECIAL_LENGTH=1000000 ;
# EXPECTED_OTHER_LENGTH=75000 ;
# LEAVE_SPECIAL=$( echo ${BIN_SIZE} ${EXPECTED_SPECIAL_LENGTH} | awk '{print $1/$2}' ) ; 
# LEAVE_OTHER=$( echo ${BIN_SIZE} ${EXPECTED_OTHER_LENGTH} | awk '{print $1/$2}' ) ;

#  for F in  "$STEP5DIR"/*bedGraph ; do
#    BASEF=$( basename $F .bedGraph ) ; 
#    echo ${BASEF} ; echo ; 
#    for VAR in LEAVE_SPECIAL LEAVE_OTHER ; do echo -e "\t${VAR}\t${!VAR}" ; done ; echo ; 
#    ${PUFFSTEP} puffcn -i ${F} --leave_special_state ${LEAVE_SPECIAL} --leave_other ${LEAVE_OTHER} --exp_decay --mu ${EMISSION_MEANS} --sigma ${EMISSION_SIGMAS} | sortBed -i - > ${BASEF}.expdecay.states.bedGraph ; 
#  done


####################################################################################################################
####################################################################################################################
## STEP 7 — COLLAPSE HMM STATE PATH BINS INTO CONTIGUOUS REGIONS WITH SAME STATE
## Take in the HMM statepath bedGraphs and collapse contiguous bins with same state into regions
## Output collapsed bedGraphs with state as the score column and length of region as an additional column
## This solves the issue of removed zero bins causing fragmentation of regions that should be contiguous.
## It can be done in the HMM step directly as well with --collapsed, but doing it as a separate step here allows for more flexibility and easier debugging/QC of the HMM state paths before collapsing.
####################################################################################################################
####################################################################################################################

mkdir -p "$STEP7DIR"
shopt -s nullglob

FILES=( "$STEP6DIR"/*.bedGraph )
(( ${#FILES[@]} )) || { echo "No files found in $STEP6DIR" >&2; exit 1; }

export COLLAPSE STEP7DIR

parallel "${PARALLEL_OPTS[@]}" '
    F={}
    BASE=$(basename "$F" .bedGraph)

    "$COLLAPSE" "$F" |
    sort -k1,1 -k2,2n --stable \
    > "$STEP7DIR/${BASE}.collapsed.bedGraph"
' ::: "${FILES[@]}"

shopt -u nullglob



## Simpler code: unparallel version that processes files sequentially but still uses collapseBedGraphRuns for the collapsing step.
# mkdir -p "$STEP7DIR"

# for F in "$STEP6DIR"/*.bedGraph ; do 
#   BASEF=$( basename ${F} .bedGraph ) ; 
#   echo ${BASEF}
#   ${COLLAPSE} ${F} | sortBed -i - > "$STEP7DIR/${BASEF}.collapsed.bedGraph"
# done



####################################################################################################################
####################################################################################################################
## STEP 8 — extract above background states.
####################################################################################################################
####################################################################################################################

mkdir -p "$STEP8DIR"
shopt -s nullglob

FILES=( "$STEP7DIR"/*.bedGraph )
(( ${#FILES[@]} )) || { echo "No files found" >&2; exit 1; }

export STEP8DIR

parallel "${PARALLEL_OPTS[@]}" '
    F={}
    B=$(basename "$F" .bedGraph)

    awk "$4>1" "$F" |
    sort -k1,1 -k2,2n --stable \
    > "$STEP8DIR/${B}.aboveBackgroundRegions.bedGraph"

    mergeBed -i "$STEP8DIR/${B}.aboveBackgroundRegions.bedGraph" \
    > "$STEP8DIR/${B}.aboveBackgroundRegions.bed"
' ::: "${FILES[@]}"

shopt -u nullglob



## simpler code: unparallel version that processes files sequentially but still uses awk and mergeBed to extract above background regions.
# mkdir -p "$STEP8DIR"

# for BDG in "$STEP7DIR"/*.bedGraph ; do 
#   BASE=$( basename ${BDG} .bedGraph ) ; 
#   echo ${BASE} ; 
#   awk '$4>1' $BDG | sortBed -i - > "$STEP8DIR/${BASE}.aboveBackgroundRegions.bedGraph" ; 
#   mergeBed -i "$STEP8DIR/${BASE}.aboveBackgroundRegions.bedGraph" > "$STEP8DIR/${BASE}.aboveBackgroundRegions.bed" ; 
# done


####################################################################################################################
####################################################################################################################
## STEP 9 — extract summit state regions.
####################################################################################################################
####################################################################################################################

mkdir -p "$STEP9DIR"
shopt -s nullglob

FILES=( "$STEP8DIR"/*.bed )
(( ${#FILES[@]} )) || { echo "ERROR: no BED files found in $STEP8DIR" >&2; exit 1; }

export STEP9DIR PUFFSTEP

parallel "${PARALLEL_OPTS[@]}" '
    BED={}
    BASE=$(basename "$BED" .bed)

    echo "Summits → $BASE" >&2

    BDG="${BED}Graph"

    if [[ ! -f "$BDG" ]]; then
        echo "ERROR: missing matching bedGraph for $BED" >&2
        exit 1
    fi

    "$PUFFSTEP" summits \
        -i "$BDG" \
        --regions "$BED" \
        --thresh_state 1 |
    sort -k1,1 -k2,2n --stable \
    > "$STEP9DIR/${BASE}.summitStateRegion.bed"
' ::: "${FILES[@]}"

shopt -u nullglob


## Simpler code: unparallel version that processes files sequentially but still uses puffStep summits to extract summit state regions from the above background regions.
# mkdir -p "$STEP9DIR"

# for BED in "$STEP8DIR"/*.bed ; do 
#   BASEF=$( basename ${BED} .bed ) ; 
#   echo ${BASEF}
#   BDG=${BED}Graph 
#   ${PUFFSTEP} summits -i ${BDG} --regions ${BED} --thresh_state 1 | sortBed -i - > "$STEP9DIR/${BASEF}.summitStateRegion.bed"
# done


####################################################################################################################
####################################################################################################################
## STEP 10 — extract summit bins.
####################################################################################################################
####################################################################################################################

mkdir -p "$STEP10DIR"
shopt -s nullglob

export STEP4DIR STEP5DIR STEP9DIR STEP10DIR PUFFSTEP PSEUDO SMOOTHING_HALFWIDTH

parallel "${PARALLEL_OPTS[@]}" '
    N={}

    echo "Processing group $N" >&2

    # --------------------------------------------------
    # SMOOTHED summit regions
    # --------------------------------------------------
    SMOOTHBDG="$STEP5DIR/${N}.TPM.medNorm.stage_${N}_median.zeroBinsRemoved.RCN_pseudo${PSEUDO}.medSmoothedHW${SMOOTHING_HALFWIDTH}.bedGraph"
    BED1=( "$STEP9DIR/${N}.TPM.medNorm.stage_${N}_median.zeroBinsRemoved.RCN_pseudo${PSEUDO}"*aboveBackgroundRegions.summitStateRegion.bed )

    if (( ${#BED1[@]} != 1 )); then
        echo "ERROR: expected 1 smoothed BED for group $N, found ${#BED1[@]}" >&2
        exit 1
    fi

    BASE=$(basename "$SMOOTHBDG" .bedGraph)

    "$PUFFSTEP" summits \
        -i "$SMOOTHBDG" \
        --regions "${BED1[0]}" |
    sort -k1,1 -k2,2n --stable \
    > "$STEP10DIR/${BASE}.summitStateRegion.bed"


    # --------------------------------------------------
    # RAW summit regions
    # --------------------------------------------------
    RAWBDG="$STEP4DIR/${N}.TPM.medNorm.stage_${N}_median.zeroBinsRemoved.RCN_pseudo${PSEUDO}.bedGraph"
    BED2=( "$STEP9DIR/${N}.TPM.medNorm.stage_${N}_median.zeroBinsRemoved.RCN_pseudo${PSEUDO}.expdecay"*".aboveBackgroundRegions.summitStateRegion.bed" )

    if (( ${#BED2[@]} != 1 )); then
        echo "ERROR: expected 1 raw BED for group $N, found ${#BED2[@]}" >&2
        exit 1
    fi

    BASE=$(basename "$RAWBDG" .bedGraph)

    "$PUFFSTEP" summits \
        -i "$RAWBDG" \
        --regions "${BED2[0]}" |
    sort -k1,1 -k2,2n --stable \
    > "$STEP10DIR/${BASE}.summitStateRegion.bed"

' ::: "${GROUP_LIST[@]:1}"

shopt -u nullglob



# ## TO OPTIONALLY ADD:
# STEP11DIR=11-medianSmoothed-log2RCN
# STEP12DIR=12-HMM-log2RCN
# STEP13DIR=13-collapsedHMM-log2RCN
# STEP14DIR=14-aboveBackground-log2RCN
# STEP15DIR=15-summitStates-log2RCN
# STEP16DIR=16-summitBins-log2RCN


## MORE CONCISE VERSION:
# parallel --keep-order --halt soon,fail=1 --jobs 0 '
# N={}
# SB="$STEP5DIR/${N}.TPM.medNorm.stage_${N}_median.zeroBinsRemoved.RCN_pseudo${PSEUDO}.medSmoothedHW${SMOOTHING_HALFWIDTH}.bedGraph"
# RB="$STEP4DIR/${N}.TPM.medNorm.stage_${N}_median.zeroBinsRemoved.RCN_pseudo${PSEUDO}.bedGraph"

# B1=( "$STEP9DIR/${N}."*summitStateRegion.bed )
# B2=( "$STEP9DIR/${N}."*expdecay*.bed )

# "$PUFFSTEP" summits -i "$SB" --regions "${B1[0]}" | sort -k1,1 -k2,2n --stable > "$STEP10DIR/$(basename "$SB" .bedGraph).summitStateRegion.bed"
# "$PUFFSTEP" summits -i "$RB" --regions "${B2[0]}" | sort -k1,1 -k2,2n --stable > "$STEP10DIR/$(basename "$RB" .bedGraph).summitStateRegion.bed"
# ' ::: "${GROUP_LIST[@]:1}"




## Simpler code: unparallel version that processes files sequentially but still uses puffStep summits to extract summit bins from the summit state regions.
# mkdir -p "$STEP10DIR"

# for N in "${GROUP_LIST[@]:1}"; do
#   BED="$STEP9DIR/${N}.TPM.medNorm.stage_${N}_median.zeroBinsRemoved.RCN_pseudo${PSEUDO}.medSmoothedHW${SMOOTHING_HALFWIDTH}*aboveBackgroundRegions.summitStateRegion.bed"
  
#   ## SMOOTHED SUMMIT BIN
#   SMOOTHBDG="$STEP5DIR/${N}.TPM.medNorm.stage_${N}_median.zeroBinsRemoved.RCN_pseudo${PSEUDO}.medSmoothedHW${SMOOTHING_HALFWIDTH}.bedGraph"
#   BASE=$( basename "${SMOOTHBDG}" .bedGraph ) ; 
#   ${PUFFSTEP} summits -i "${SMOOTHBDG}" --regions "${BED}" > "${STEP10DIR}/${BASE}.summitStateRegion.bed"

#   ## UNSMOOTHED SUMMIT BIN
#   RAWBDG="$STEP4DIR/${N}.TPM.medNorm.stage_${N}_median.zeroBinsRemoved.RCN_pseudo${PSEUDO}.bedGraph"
#   BASE=$( basename "${RAWBDG}" .bedGraph ) ; 
#   BED="$STEP9DIR/${N}.TPM.medNorm.stage_${N}_median.zeroBinsRemoved.RCN_pseudo${PSEUDO}.expdecay*.aboveBackgroundRegions.summitStateRegion.bed"
#   ${PUFFSTEP} summits -i "${RAWBDG}" --regions "${BED}" > "${STEP10DIR}/${BASE}.summitStateRegion.bed"
# done





