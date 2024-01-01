# Nextflow command to pre-process scRNA-seq raw data

# scRNA-seq raw data: fastq.gz files
# Pipeline: nf-core/scrnaseq (https://nf-co.re/scrnaseq), Version: 2.3.0
# Input csv file format: sample[sample_name], fastq_1[fastq.gz], fastq_2[fastq.gz]
# Reference genome: GRCh38 Gencode v32, download command (provided by Cellranger): wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

nextflow run nf-core/scrnaseq --input sample_scrnaseq_config.csv -profile docker \
  --aligner cellranger -r 2.3.0 --outdir output_cellranger --igenomes_ignore \
  --fasta refdata-gex-GRCh38-2020-A/fasta/genome.fa --gtf refdata-gex-GRCh38-2020-A/genes/genes.gtf
  