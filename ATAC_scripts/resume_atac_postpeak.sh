#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/d/lxk/project/lishan-20260613"
ATAC="$PROJECT/analysis/ATAC"
SAMPLES="$ATAC/metadata/samples.tsv"
REF="$HOME/database/source/human"
GENOME="$REF/genome.fa"
GTF="$REF/genes.filtered.gtf"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate atac_py

mkdir -p "$ATAC"/{logs,qc,counts,annotation,figures,matrix}
exec > >(tee -a "$ATAC/logs/resume_postpeak.$(date +%Y%m%d_%H%M%S).log") 2>&1
echo "[RESUME] $(date)"

if [ ! -s "$GENOME.fai" ]; then
  samtools faidx "$GENOME"
fi
awk '{print $1"\t"$2}' "$GENOME.fai" > "$REF/genome.chrom.sizes"

python "$ATAC/scripts/prepare_scrna_gene_sets.py"
python "$ATAC/scripts/make_tss_bed.py" "$GTF" "$ATAC/annotation/genes.tss.unsorted.bed" "$ATAC/annotation/genes.body.unsorted.bed"
bedtools sort -faidx "$GENOME.fai" -i "$ATAC/annotation/genes.tss.unsorted.bed" > "$ATAC/annotation/genes.tss.bed"
bedtools sort -faidx "$GENOME.fai" -i "$ATAC/annotation/genes.body.unsorted.bed" > "$ATAC/annotation/genes.body.bed"

echo "[CONSENSUS] merge narrowPeaks"
cat "$ATAC"/peaks/*/*_peaks.narrowPeak \
  | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ {print $1,$2,$3}' \
  | bedtools sort -faidx "$GENOME.fai" -i - \
  | bedtools merge -i - \
  > "$ATAC/peaks/consensus_peaks.bed"

echo "[COUNTS] bedtools multicov"
mapfile -t bam_files < <(tail -n +2 "$SAMPLES" | awk -v b="$ATAC/bam" '{print b"/"$1".final.bam"}')
bedtools multicov -bams "${bam_files[@]}" -bed "$ATAC/peaks/consensus_peaks.bed" > "$ATAC/counts/consensus_peak_counts.raw.tsv"
{
  printf "chrom\tstart\tend"
  tail -n +2 "$SAMPLES" | awk '{printf "\t"$1}'
  printf "\n"
  cat "$ATAC/counts/consensus_peak_counts.raw.tsv"
} > "$ATAC/counts/consensus_peak_counts.tsv"

echo "[ANNOTATION] nearest genes and scRNA links"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"peak_"NR}' "$ATAC/peaks/consensus_peaks.bed" > "$ATAC/peaks/consensus_peaks.named.bed"
bedtools closest -sorted -g "$REF/genome.chrom.sizes" -d -a "$ATAC/peaks/consensus_peaks.named.bed" -b "$ATAC/annotation/genes.body.bed" \
  > "$ATAC/annotation/consensus_peaks_nearest_gene.raw.tsv"
python "$ATAC/scripts/annotate_and_link.py"

echo "[QC] sample metrics"
{
  printf "sample\tgroup\tfinal_mapped_fragments\tpeaks\tfrip\n"
  tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r sample group r1 r2; do
    bam="$ATAC/bam/${sample}.final.bam"
    reads=$(samtools view -@ 4 -c "$bam")
    fragments=$((reads / 2))
    peaks=$(wc -l < "$ATAC/peaks/${sample}/${sample}_peaks.narrowPeak")
    in_peak_reads=$(bedtools intersect -u -a "$bam" -b "$ATAC/peaks/consensus_peaks.bed" | samtools view -c -)
    frip=$(awk -v a="$in_peak_reads" -v b="$reads" 'BEGIN{if(b>0) printf "%.5f", a/b; else print "NA"}')
    printf "%s\t%s\t%s\t%s\t%s\n" "$sample" "$group" "$fragments" "$peaks" "$frip"
  done
} > "$ATAC/qc/sample_qc_metrics.tsv"

echo "[DEEPTOOLS] TSS matrix/profile"
computeMatrix reference-point \
  --referencePoint TSS \
  -b 3000 -a 3000 \
  -R "$ATAC/annotation/genes.tss.bed" \
  -S "$ATAC"/bigwig/NP1.CPM.bw "$ATAC"/bigwig/NP2.CPM.bw "$ATAC"/bigwig/NP3.CPM.bw "$ATAC"/bigwig/PE1.CPM.bw \
  --samplesLabel NP1 NP2 NP3 PE1 \
  --skipZeros \
  -p 8 \
  -o "$ATAC/matrix/tss_CPM_matrix.gz"
plotProfile -m "$ATAC/matrix/tss_CPM_matrix.gz" -out "$ATAC/figures/ATAC_TSS_profile.pdf" --colors "#4C78A8" "#4C78A8" "#4C78A8" "#D55E00"
plotHeatmap -m "$ATAC/matrix/tss_CPM_matrix.gz" -out "$ATAC/figures/ATAC_TSS_heatmap.pdf" --colorMap RdBu_r --whatToShow "heatmap and colorbar"

echo "[FIGURE] summary PDF"
python "$ATAC/scripts/plot_atac_summary.py"
echo "[DONE] $(date)"
