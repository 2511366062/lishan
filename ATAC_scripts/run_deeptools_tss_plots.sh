#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/d/lxk/project/lishan-20260613"
ATAC="$PROJECT/analysis/ATAC"
OUT="$PROJECT/fig/ATAC_deepTools"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate atac_py

mkdir -p "$ATAC/matrix" "$OUT"

awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ && $4 !~ /^ENSG/ {print}' "$ATAC/annotation/genes.tss.bed" \
  | awk 'NR % 3 == 1 {print}' \
  | head -n 12000 \
  | sort -k1,1 -k2,2n \
  > "$ATAC/annotation/genes.tss.deepTools_12k.bed"

computeMatrix reference-point \
  --referencePoint TSS \
  -b 3000 -a 3000 \
  -R "$ATAC/annotation/genes.tss.deepTools_12k.bed" \
  -S "$ATAC/bigwig/NP1.CPM.bw" "$ATAC/bigwig/NP2.CPM.bw" "$ATAC/bigwig/NP3.CPM.bw" "$ATAC/bigwig/PE1.CPM.bw" \
  --samplesLabel NP1 NP2 NP3 PE1 \
  --skipZeros \
  --missingDataAsZero \
  -p 8 \
  -o "$ATAC/matrix/deepTools_TSS_12k_CPM_matrix.gz" \
  --outFileSortedRegions "$ATAC/matrix/deepTools_TSS_12k_sorted_regions.bed"

plotHeatmap \
  -m "$ATAC/matrix/deepTools_TSS_12k_CPM_matrix.gz" \
  -out "$OUT/deepTools_TSS_heatmap_profile.pdf" \
  --plotTitle "ATAC-seq signal around TSS" \
  --refPointLabel "TSS" \
  --regionsLabel "genes" \
  --samplesLabel NP1 NP2 NP3 PE1 \
  --colorMap RdYlBu \
  --whatToShow "plot, heatmap and colorbar" \
  --heatmapHeight 12 \
  --heatmapWidth 3.0 \
  --xAxisLabel "" \
  --startLabel "-3 kb" \
  --endLabel "3 kb" \
  --sortUsing mean \
  --sortRegions descend \
  --zMin 0 \
  --zMax 1.2

plotProfile \
  -m "$ATAC/matrix/deepTools_TSS_12k_CPM_matrix.gz" \
  -out "$OUT/deepTools_TSS_profile.pdf" \
  --plotTitle "Average ATAC-seq signal around TSS" \
  --refPointLabel "TSS" \
  --regionsLabel "genes" \
  --samplesLabel NP1 NP2 NP3 PE1 \
  --colors "#4C78A8" "#6BAED6" "#9ECAE1" "#D55E00" \
  --plotHeight 4 \
  --plotWidth 5

echo "Wrote deepTools plots to $OUT"
