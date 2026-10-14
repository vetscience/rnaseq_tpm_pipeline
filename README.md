# rnaseq_tpm_pipeline
Python pipeline for gene-level RNA-seq quantification across multiple species. Converts SRA runs to FASTQ, checks quality (FastQC), trims paired-end reads (Trimmomatic), aligns with HISAT2 or, for long reads, minimap2, counts reads per gene (featureCounts) and exports per-sample TPM tables as CSV/Excel. Commands and tool versions are logged.
