#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/d/lxk/project/lishan-20260613"
ATAC="$PROJECT/analysis/ATAC"
SAMPLES="$ATAC/metadata/samples.tsv"
REF="$HOME/database/source/human"
GENOME="$REF/genome.fa"
GTF="$REF/genes.filtered.gtf"
INDEX_PREFIX="$REF/bowtie2/genome"
THREADS="${THREADS:-20}"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate atac_py

mkdir -p "$ATAC"/{logs,qc,trimmed,bam,peaks,bigwig,counts,annotation,figures,matrix} "$REF/bowtie2"

exec > >(tee -a "$ATAC/logs/run_bulk_atac.$(date +%Y%m%d_%H%M%S).log") 2>&1
echo "[START] $(date)"
echo "[ENV] $(which python)"
echo "[THREADS] $THREADS"

python "$ATAC/scripts/prepare_scrna_gene_sets.py"

if [ ! -s "$GENOME.fai" ]; then
  echo "[REF] samtools faidx"
  samtools faidx "$GENOME"
fi

awk '{print $1"\t"$2}' "$GENOME.fai" > "$REF/genome.chrom.sizes"
awk '$1 ~ /^chr([0-9]+|X|Y)$/ {print $1"\t0\t"$2}' "$REF/genome.chrom.sizes" > "$REF/genome.primary.bed"

if [ ! -s "${INDEX_PREFIX}.1.bt2" ] && [ ! -s "${INDEX_PREFIX}.1.bt2l" ]; then
  echo "[REF] bowtie2-build"
  bowtie2-build --threads "$THREADS" "$GENOME" "$INDEX_PREFIX"
fi

python "$ATAC/scripts/make_tss_bed.py" "$GTF" "$ATAC/annotation/genes.tss.unsorted.bed" "$ATAC/annotation/genes.body.unsorted.bed"
bedtools sort -faidx "$GENOME.fai" -i "$ATAC/annotation/genes.tss.unsorted.bed" > "$ATAC/annotation/genes.tss.bed"
bedtools sort -faidx "$GENOME.fai" -i "$ATAC/annotation/genes.body.unsorted.bed" > "$ATAC/annotation/genes.body.bed"

tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r sample group r1 r2; do
  echo "[SAMPLE] $sample $group"
  r1_abs="$PROJECT/$r1"
  r2_abs="$PROJECT/$r2"
  trim_r1="$ATAC/trimmed/${sample}_R1.fastp.fq.gz"
  trim_r2="$ATAC/trimmed/${sample}_R2.fastp.fq.gz"
  raw_bam="$ATAC/bam/${sample}.raw.bam"
  name_bam="$ATAC/bam/${sample}.name.bam"
  fix_bam="$ATAC/bam/${sample}.fixmate.bam"
  sort_bam="$ATAC/bam/${sample}.sorted.bam"
  dedup_bam="$ATAC/bam/${sample}.dedup.bam"
  final_bam="$ATAC/bam/${sample}.final.bam"

  if [ ! -s "$trim_r1" ]; then
    fastp \
      -i "$r1_abs" -I "$r2_abs" \
      -o "$trim_r1" -O "$trim_r2" \
      --detect_adapter_for_pe \
      --thread 8 \
      --html "$ATAC/qc/${sample}.fastp.html" \
      --json "$ATAC/qc/${sample}.fastp.json"
  fi

  if [ ! -s "$final_bam" ]; then
    bowtie2 --very-sensitive -X 2000 -p "$THREADS" -x "$INDEX_PREFIX" -1 "$trim_r1" -2 "$trim_r2" \
      2> "$ATAC/logs/${sample}.bowtie2.log" \
      | samtools view -@ 4 -b -q 30 -f 2 -F 1804 -o "$raw_bam" -
    samtools sort -@ "$THREADS" -n -o "$name_bam" "$raw_bam"
    samtools fixmate -@ "$THREADS" -m "$name_bam" "$fix_bam"
    samtools sort -@ "$THREADS" -o "$sort_bam" "$fix_bam"
    samtools markdup -@ "$THREADS" -r "$sort_bam" "$dedup_bam"
    samtools view -@ 4 -h "$dedup_bam" \
      | awk 'BEGIN{OFS="\t"} /^@/ || ($3!="chrM" && $3!="MT" && $3!="M")' \
      | samtools view -@ 4 -b -o "$final_bam" -
    samtools index -@ "$THREADS" "$final_bam"
    rm -f "$raw_bam" "$name_bam" "$fix_bam" "$sort_bam" "$dedup_bam"
  fi

  samtools flagstat -@ "$THREADS" "$final_bam" > "$ATAC/qc/${sample}.flagstat.txt"
  samtools idxstats "$final_bam" > "$ATAC/qc/${sample}.idxstats.txt"

  if [ ! -s "$ATAC/peaks/${sample}/${sample}_peaks.narrowPeak" ]; then
    mkdir -p "$ATAC/peaks/${sample}"
    macs3 callpeak \
      -t "$final_bam" \
      -f BAMPE \
      -g hs \
      -n "$sample" \
      --outdir "$ATAC/peaks/${sample}" \
      -q 0.01 \
      --keep-dup all \
      --call-summits
  fi

  if [ ! -s "$ATAC/bigwig/${sample}.CPM.bw" ]; then
    bamCoverage \
      -b "$final_bam" \
      -o "$ATAC/bigwig/${sample}.CPM.bw" \
      --normalizeUsing CPM \
      --binSize 10 \
      --extendReads \
      -p 8
  fi
done

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
if [ ! -s "$ATAC/matrix/tss_CPM_matrix.gz" ]; then
  computeMatrix reference-point \
    --referencePoint TSS \
    -b 3000 -a 3000 \
    -R "$ATAC/annotation/genes.tss.bed" \
    -S "$ATAC"/bigwig/*.CPM.bw \
    --skipZeros \
    -p 8 \
    -o "$ATAC/matrix/tss_CPM_matrix.gz"
fi
plotProfile -m "$ATAC/matrix/tss_CPM_matrix.gz" -out "$ATAC/figures/ATAC_TSS_profile.pdf" --perGroup --colors "#4C78A8" "#4C78A8" "#4C78A8" "#D55E00"
plotHeatmap -m "$ATAC/matrix/tss_CPM_matrix.gz" -out "$ATAC/figures/ATAC_TSS_heatmap.pdf" --colorMap RdBu_r --whatToShow "heatmap and colorbar"

echo "[FIGURE] summary PDF"
python "$ATAC/scripts/plot_atac_summary.py"

echo "[DONE] $(date)"
