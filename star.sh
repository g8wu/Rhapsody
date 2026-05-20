#!/bin/bash

files="G:/bulkRNAseq/RAWdatabulkRNAseq"
genome="G:/bulkRNAseq/reference_files"
gtf="G:/bulkRNAseq/reference_files/Homo_sapiens.GRCh38.99.gtf"
rseqcBed="G:/bulkRNAseq/reference_files"
outDir="G:/bulkRNAseq/star_out"

#mkdir -p "$outDir"/{star,counts,rseqc}

# Run STAR in read pairs
find "$files" -type f -name "*_1.fq.gz" | while read r1; do
  # Get file base name
  sample=$(basename "$r1" _1.fastq.gz)
  r2="${r1/_1/_2}"
  
  echo "Processing $sample ________________________________"
  # Run alignment
  STAR --runThreadN 16 \
       --genomeDir "$genome" \
       --readFilesIn "$r1" "$r2" \
       --readFilesCommand zcat \
       --outFileNamePrefix "$outDir/star/${sample}_" \
       --outSAMtype BAM SortedByCoordinate

  BAM="$outDir/star/${sample}_.bam"

  # featureCounts quant
  echo "Starting FeatureCounts: $sample ........."

  featureCounts -T 16 -p -t exon -g gene_id \
  -a "$gtf" \
  -o "$outDir/counts/${sample}_counts.txt" \
  "$BAM"

  echo "Starting RseQC: $sample ........."
  
  # RSeqc: geneBody_coverage
  geneBody_coverage.py -r "$rseqcBed" -i "$BAM" -o "$outDir/rseqc/${sample}_geneBody"

  # RSeqc: read_distribution
  read_distribution.py -r "$rseqcBed" -i "$BAM" > "$outDir/rseqc/${sample}_readDistr.txt"
done
