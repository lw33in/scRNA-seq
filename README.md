# scRNA-seq
- **single-cell RNA-sequencing (scRNA-seq)** is a genomic approach for the detection and quantitative analysis of mRNA molecules in a biological sample and is useful for studying cellular responses.
- **Input cell requirement**:
  * ~1.2M cells per tube with > 80% cell viability
- **Data type**:
  * Paired-end 150bp
  * ~60G/sample
- **Library construction**:
  * [Chromium Single Cell 3' Reagent Kits](https://www.10xgenomics.com/support/single-cell-gene-expression/documentation/steps/library-prep/chromium-single-cell-3-reagent-kits-user-guide-v-3-1-chemistry-dual-index)
- **Sequencing platform**:
  * Illumina® NovaSeq 6000


## scRNA-seq Experimental and Analysis Workflow

> Note: This is a standard workflow designed for scRNA-seq. Based on the specific biological questions and datasets, modifications of this workflow might be needed.

- **Pre-processing**:
    * [Cellranger from 10X Genomics](https://www.10xgenomics.com/support/software/cell-ranger) is used for conducting raw data demultiplexing, barcoding process, gene counting, transcript assembly, annotation, and barcode analysis.
    * In this workflow, [nf-core/scrnaseq/2.3.0](https://nf-co.re/scrnaseq/2.3.0) is used to run Cellranger. v2.3.0 is an older version of the pipeline. If you wish to use the latest release, please visit https://nf-co.re/scrnaseq/
- **Downstream analysis**:
  * [Seurat](https://satijalab.org/seurat/) is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. Seurat is developed and maintained by the [Satija lab](https://satijalab.org/seurat/) and is released under the MIT license.

<p align="center">
<img width="463" alt="image" src="https://github.com/stephniw/scRNA-seq/assets/120678930/663e684b-4f35-49df-8ac1-1288b20f88d7">
</p>

- **Reference Genome**:
    * Gencode v32 reference genome
    * Command line download: <code>wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz</code>
    * Direct download: [Human reference (GRCh38) dataset required for Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest/)

## RNA-seq Computational Analysis
- `01_nextflow_scRNAseq.sh` General command lines used to process raw reads using nf-core pipelines.
  * sample_scrnaseq_config.csv:
  ```
  sample,fastq_1,fastq_2,expected_cells
  ctrl,ctrl_S1_L007_R1_001.fastq.gz,ctrl_S1_L007_R2_001.fastq.gz,"10000"
  exp,exp_S1_L008_R1_001.fastq.gz,exp_S1_L008_R2_001.fastq.gz,"10000"
  ```

  ```
  nextflow run nf-core/scrnaseq --input sample_scrnaseq_config.csv -profile docker \
  --aligner cellranger -r 2.3.0 --outdir output_cellranger --igenomes_ignore \
  --fasta refdata-gex-GRCh38-2020-A/fasta/genome.fa --gtf refdata-gex-GRCh38-2020-A/genes/genes.gtf
  ```
    > Note If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with -profile test before running the workflow on actual data.
  
- `02_QC.R` Perform QC on raw count data, filter out poor quality cells and unwanted variances, generate associated plots.
- `03_integration.R` Normalize, transform, and integrate datasets across experimental batches or conditions.
- `04_clustering.R` Perform PCA and unsupervised clustering to assign cells to clusters.
- `05_annotation.R` Compare expression profiles of clusters to previously annotated reference datasets.
- `06_pseudobulk_DEA.R` Perform pseudobulk differential expression analysis across different clusters / cell types.
- `07_cluster_annotation.R` Summarize cell origins in different clusters / cell types.

## Tools 
- R-4.3.2
- nextflow 2.3.0
- Seurat v5
  
