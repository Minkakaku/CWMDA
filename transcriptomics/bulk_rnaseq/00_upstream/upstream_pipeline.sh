#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  upstream_pipeline.sh [options]

Options:
  --sample-sheet PATH     Sample sheet with columns: sample_id, layout, sra (default: data/transcriptome/metadata/sample_sheet.tsv)
  --acc-list PATH         Plain text list of SRA accessions (used when sample sheet is unavailable)
  --default-layout LAYOUT Default layout when using --acc-list (SE or PE, default: SE)
  --aligner NAME          Aligner: star or hisat2 (default: star)
  --quant NAME            Quantifier: featurecounts or rsem (default: featurecounts)
  --genome-dir PATH       STAR genomeDir (required for STAR)
  --hisat2-index PATH     HISAT2 index base path (required for HISAT2)
  --gtf PATH              GTF annotation (required for featureCounts)
  --rsem-ref PATH         RSEM reference prefix (required for RSEM)
  --threads N             Threads (default: 8)
  --work-dir PATH         Working directory for outputs (default: current directory)
  --skip-download         Skip SRA prefetch/fastq-dump
  --skip-qc               Skip fastp QC/trim
  --help                  Show this help

Examples:
  ./upstream_pipeline.sh --aligner star --quant featurecounts \
    --genome-dir /data/index/STAR --gtf /data/genes.gtf

  ./upstream_pipeline.sh --aligner hisat2 --quant rsem \
    --hisat2-index /data/index/hisat2/genome --rsem-ref /data/index/rsem/ref
USAGE
}

log() {
  printf "[%s] %s\n" "$(date '+%F %T')" "$*"
}

die() {
  printf "[ERROR] %s\n" "$*" >&2
  exit 1
}

command_exists() {
  command -v "$1" >/dev/null 2>&1
}

ALIGNER="star"
QUANT="featurecounts"
THREADS=8
WORK_DIR="$(pwd)"
SKIP_DOWNLOAD=false
SKIP_QC=false
DEFAULT_LAYOUT="SE"

BASE_DIR=$(git -C "$WORK_DIR" rev-parse --show-toplevel 2>/dev/null || true)
SAMPLE_SHEET_DEFAULT="${BASE_DIR:-$WORK_DIR}/data/transcriptome/metadata/sample_sheet.tsv"
SAMPLE_SHEET="$SAMPLE_SHEET_DEFAULT"
ACC_LIST=""
GENOME_DIR=""
HISAT2_INDEX=""
GTF=""
RSEM_REF=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-sheet)
      SAMPLE_SHEET="$2"; shift 2;;
    --acc-list)
      ACC_LIST="$2"; shift 2;;
    --default-layout)
      DEFAULT_LAYOUT="$2"; shift 2;;
    --aligner)
      ALIGNER="$2"; shift 2;;
    --quant)
      QUANT="$2"; shift 2;;
    --genome-dir)
      GENOME_DIR="$2"; shift 2;;
    --hisat2-index)
      HISAT2_INDEX="$2"; shift 2;;
    --gtf)
      GTF="$2"; shift 2;;
    --rsem-ref)
      RSEM_REF="$2"; shift 2;;
    --threads)
      THREADS="$2"; shift 2;;
    --work-dir)
      WORK_DIR="$2"; shift 2;;
    --skip-download)
      SKIP_DOWNLOAD=true; shift;;
    --skip-qc)
      SKIP_QC=true; shift;;
    --help)
      usage; exit 0;;
    *)
      die "Unknown argument: $1";;
  esac
 done

case "${ALIGNER,,}" in
  star|hisat2) ;; 
  *) die "--aligner must be star or hisat2";;
 esac

case "${QUANT,,}" in
  featurecounts|rsem) ;;
  *) die "--quant must be featurecounts or rsem";;
 esac

if [[ -z "$ACC_LIST" && ! -f "$SAMPLE_SHEET" ]]; then
  die "Sample sheet not found: $SAMPLE_SHEET (or provide --acc-list)"
fi

if [[ -z "$ACC_LIST" ]]; then
  if ! awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="layout") l=1; if (!l) exit 1}' "$SAMPLE_SHEET"; then
    die "Sample sheet missing 'layout' column: $SAMPLE_SHEET"
  fi
fi

if [[ "${ALIGNER,,}" == "star" && -z "$GENOME_DIR" ]]; then
  die "STAR requires --genome-dir"
fi

if [[ "${ALIGNER,,}" == "hisat2" && -z "$HISAT2_INDEX" ]]; then
  die "HISAT2 requires --hisat2-index"
fi

if [[ "${QUANT,,}" == "featurecounts" && -z "$GTF" ]]; then
  die "featureCounts requires --gtf"
fi

if [[ "${QUANT,,}" == "rsem" && -z "$RSEM_REF" ]]; then
  die "RSEM requires --rsem-ref"
fi

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

RAW_DIR="00.rawdata"
FASTQ_DIR="01.fastq"
QC_DIR="02.cleandata"
ALIGN_DIR="03.align"
QUANT_DIR="04.quant"
LOG_DIR="logs"

mkdir -p "$RAW_DIR" "$FASTQ_DIR" "$QC_DIR" "$ALIGN_DIR" "$QUANT_DIR" "$LOG_DIR"

get_sra_by_layout() {
  local target_layout=$1
  awk -F'\t' -v target="$target_layout" '
    NR==1 {
      for (i=1; i<=NF; i++) {
        if ($i=="layout") layout=i
        if ($i=="sra") sra=i
      }
      next
    }
    toupper($layout)==toupper(target) {print $sra}
  ' "$SAMPLE_SHEET"
}

get_all_sra() {
  if [[ -n "$ACC_LIST" ]]; then
    cat "$ACC_LIST"
    return 0
  fi
  awk -F'\t' '
    NR==1 {
      for (i=1; i<=NF; i++) {
        if ($i=="sra") sra=i
      }
      next
    }
    {print $sra}
  ' "$SAMPLE_SHEET"
}

collect_layouts() {
  if [[ -n "$ACC_LIST" ]]; then
    case "${DEFAULT_LAYOUT^^}" in
      SE)
        [[ "${1^^}" == "SE" ]] && get_all_sra || true
        ;;
      PE)
        [[ "${1^^}" == "PE" ]] && get_all_sra || true
        ;;
      *) die "--default-layout must be SE or PE";;
    esac
    return 0
  fi
  get_sra_by_layout "$1"
}

require_tools=(prefetch fastq-dump fastp samtools)
if [[ "${ALIGNER,,}" == "star" ]]; then
  require_tools+=(STAR)
else
  require_tools+=(hisat2)
fi
if [[ "${QUANT,,}" == "featurecounts" ]]; then
  require_tools+=(featureCounts)
else
  require_tools+=(rsem-calculate-expression)
fi

for tool in "${require_tools[@]}"; do
  if ! command_exists "$tool"; then
    die "Missing required tool: $tool"
  fi
 done

if [[ "$SKIP_DOWNLOAD" == false ]]; then
  log "Downloading SRA data"
  ACC_LIST_PATH="$ACC_LIST"
  if [[ -z "$ACC_LIST_PATH" ]]; then
    ACC_LIST_PATH="$WORK_DIR/acc_list.txt"
    get_all_sra > "$ACC_LIST_PATH"
  fi
  prefetch --option-file "$ACC_LIST_PATH" -O "$RAW_DIR" 2> "$LOG_DIR/prefetch.err" > "$LOG_DIR/prefetch.log"

  log "Converting SRA to FASTQ"
  while read -r id; do
    fastq-dump "$RAW_DIR/$id/$id.sra" --gzip --split-e --defline-seq '@$ac-$si/$ri' --defline-qual '+' \
      -O "$FASTQ_DIR" 2> "$LOG_DIR/${id}.fastq.err" > "$LOG_DIR/${id}.fastq.log"
  done < "$ACC_LIST_PATH"
else
  log "Skipping download/fastq conversion"
fi

if [[ "$SKIP_QC" == false ]]; then
  log "Running fastp QC"
  while read -r id; do
    fastp -w "$THREADS" -i "$FASTQ_DIR/${id}.fastq.gz" -o "$QC_DIR/${id}.clean.fq.gz" \
      -j "$QC_DIR/${id}.fastp.json" -h "$QC_DIR/${id}.fastp.html" \
      2> "$LOG_DIR/${id}.fastp.err" > "$LOG_DIR/${id}.fastp.log"
  done < <(collect_layouts "SE")

  while read -r id; do
    fastp -w "$THREADS" -i "$FASTQ_DIR/${id}_1.fastq.gz" -I "$FASTQ_DIR/${id}_2.fastq.gz" \
      -o "$QC_DIR/${id}_1.clean.fq.gz" -O "$QC_DIR/${id}_2.clean.fq.gz" \
      -j "$QC_DIR/${id}.fastp.json" -h "$QC_DIR/${id}.fastp.html" \
      2> "$LOG_DIR/${id}.fastp.err" > "$LOG_DIR/${id}.fastp.log"
  done < <(collect_layouts "PE")
else
  log "Skipping QC"
fi

log "Aligning reads with ${ALIGNER^^}"
if [[ "${ALIGNER,,}" == "star" ]]; then
  STAR_COMMON=(--runThreadN "$THREADS" --readFilesCommand zcat --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --genomeDir "$GENOME_DIR")
  while read -r id; do
    STAR "${STAR_COMMON[@]}" --readFilesIn "$QC_DIR/${id}.clean.fq.gz" \
      --outFileNamePrefix "$ALIGN_DIR/${id}." 2> "$LOG_DIR/${id}.star.err" > "$LOG_DIR/${id}.star.log"
    samtools index -@ "$THREADS" "$ALIGN_DIR/${id}.Aligned.sortedByCoord.out.bam" \
      2> "$LOG_DIR/${id}.index.err" > "$LOG_DIR/${id}.index.log"
  done < <(collect_layouts "SE")

  while read -r id; do
    STAR "${STAR_COMMON[@]}" --readFilesIn "$QC_DIR/${id}_1.clean.fq.gz" "$QC_DIR/${id}_2.clean.fq.gz" \
      --outFileNamePrefix "$ALIGN_DIR/${id}." 2> "$LOG_DIR/${id}.star.err" > "$LOG_DIR/${id}.star.log"
    samtools index -@ "$THREADS" "$ALIGN_DIR/${id}.Aligned.sortedByCoord.out.bam" \
      2> "$LOG_DIR/${id}.index.err" > "$LOG_DIR/${id}.index.log"
  done < <(collect_layouts "PE")
else
  while read -r id; do
    hisat2 -p "$THREADS" -x "$HISAT2_INDEX" -U "$QC_DIR/${id}.clean.fq.gz" \
      2> "$LOG_DIR/${id}.hisat2.err" | \
      samtools sort -@ "$THREADS" -o "$ALIGN_DIR/${id}.sorted.bam" - \
      2> "$LOG_DIR/${id}.sort.err" > "$LOG_DIR/${id}.sort.log"
    samtools index -@ "$THREADS" "$ALIGN_DIR/${id}.sorted.bam" \
      2> "$LOG_DIR/${id}.index.err" > "$LOG_DIR/${id}.index.log"
  done < <(collect_layouts "SE")

  while read -r id; do
    hisat2 -p "$THREADS" -x "$HISAT2_INDEX" -1 "$QC_DIR/${id}_1.clean.fq.gz" -2 "$QC_DIR/${id}_2.clean.fq.gz" \
      2> "$LOG_DIR/${id}.hisat2.err" | \
      samtools sort -@ "$THREADS" -o "$ALIGN_DIR/${id}.sorted.bam" - \
      2> "$LOG_DIR/${id}.sort.err" > "$LOG_DIR/${id}.sort.log"
    samtools index -@ "$THREADS" "$ALIGN_DIR/${id}.sorted.bam" \
      2> "$LOG_DIR/${id}.index.err" > "$LOG_DIR/${id}.index.log"
  done < <(collect_layouts "PE")
fi

log "Quantifying with ${QUANT^^}"
if [[ "${QUANT,,}" == "featurecounts" ]]; then
  while read -r id; do
    bam_path="$ALIGN_DIR/${id}.Aligned.sortedByCoord.out.bam"
    if [[ ! -f "$bam_path" ]]; then
      bam_path="$ALIGN_DIR/${id}.sorted.bam"
    fi
    featureCounts -T "$THREADS" -t exon -g gene_id -a "$GTF" -o "$QUANT_DIR/${id}.featureCounts.txt" "$bam_path" \
      2> "$LOG_DIR/${id}.featureCounts.err" > "$LOG_DIR/${id}.featureCounts.log"
  done < <(collect_layouts "SE")

  while read -r id; do
    bam_path="$ALIGN_DIR/${id}.Aligned.sortedByCoord.out.bam"
    if [[ ! -f "$bam_path" ]]; then
      bam_path="$ALIGN_DIR/${id}.sorted.bam"
    fi
    featureCounts -T "$THREADS" -p -t exon -g gene_id -a "$GTF" -o "$QUANT_DIR/${id}.featureCounts.txt" "$bam_path" \
      2> "$LOG_DIR/${id}.featureCounts.err" > "$LOG_DIR/${id}.featureCounts.log"
  done < <(collect_layouts "PE")
else
  while read -r id; do
    bam_path="$ALIGN_DIR/${id}.Aligned.sortedByCoord.out.bam"
    if [[ ! -f "$bam_path" ]]; then
      bam_path="$ALIGN_DIR/${id}.sorted.bam"
    fi
    rsem-calculate-expression --bam -p "$THREADS" "$bam_path" "$RSEM_REF" "$QUANT_DIR/${id}" \
      2> "$LOG_DIR/${id}.rsem.err" > "$LOG_DIR/${id}.rsem.log"
  done < <(collect_layouts "SE")

  while read -r id; do
    bam_path="$ALIGN_DIR/${id}.Aligned.sortedByCoord.out.bam"
    if [[ ! -f "$bam_path" ]]; then
      bam_path="$ALIGN_DIR/${id}.sorted.bam"
    fi
    rsem-calculate-expression --bam --paired-end -p "$THREADS" "$bam_path" "$RSEM_REF" "$QUANT_DIR/${id}" \
      2> "$LOG_DIR/${id}.rsem.err" > "$LOG_DIR/${id}.rsem.log"
  done < <(collect_layouts "PE")
fi

log "Pipeline finished"
