# rnaseq_tpm_pipeline
Python pipeline for gene-level RNA-seq quantification across multiple species. Converts SRA runs to FASTQ, checks quality (FastQC), trims paired-end reads (Trimmomatic), aligns with HISAT2 or, for long reads, minimap2, counts reads per gene (featureCounts) and exports per-sample TPM tables as CSV/Excel. Commands and tool versions are logged.

This code accompanies the following article. If you use it, please cite:
Kaethner, M., Young, N. D., Sleebs, B. E., Lee, E. F., Fairlie, W. D., & Gasser, R. B. Comparative reconstruction of intrinsic apoptosis in clade V nematodes reveals a conserved execution module and heterogeneous regulatory architecture. *Computational and Structural Biotechnology Journal* (2026). https://doi.org/10.34133/csbj.0239

## Workflow (per species)

| Step | Tool |
|---|---|
| 1. Reference preparation: unpack genome/annotation (`.zip`/`.gz`), convert GFF3 → GTF | gffread |
| 2. Genome index (short-read species only) | hisat2-build |
| 3. SRA → FASTQ | fasterq-dump (SRA Toolkit) |
| 4. Raw-read quality control | FastQC |
| 5. Adapter and quality trimming (paired-end) | Trimmomatic |
| 6. Alignment, sorting, indexing | HISAT2 (short reads) or minimap2 `-ax splice:hq` (long reads); samtools |
| 7. Trimmed-read quality control | FastQC |
| 8. Merge all runs of a species (`Merged_Total` column) | samtools merge |
| 9. Gene-level counting (exons, `gene_id`, unstranded) | featureCounts (Subread) |
| 10. TPM normalisation | Python / pandas |

TPM is calculated as

```
RPK_i = count_i / (length_i / 1000)
TPM_i = RPK_i / sum_j(RPK_j) * 1,000,000
```

where `length` is the featureCounts `Length` column (non-overlapping exonic bases of the gene).

## Installation

Requires Linux and conda (e.g. [Miniforge](https://github.com/conda-forge/miniforge)).

Download `rnaseq_tpm_pipeline.py` and `environment.yml` from this repository.
In the folder containing them, create and activate the environment:

```bash
conda env create -f environment.yml
conda activate rnaseq_tpm_env
```

Copy the Trimmomatic adapter file (installed with the environment) next to the script:

```bash
cp $CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa .
```

## Input

All paths are relative to the script folder:

```
rnaseq_tpm_pipeline.py
TruSeq3-PE.fa
SRR_to_run/
    Necator_americanus.txt
    ...
transcriptomes/
    Necator_americanus/
        <genome>        *genomic.fna | *genomic.fa | *_genome           (may be .gz)
        <annotation>    *.gtf | *.gff3 | *.gff | *_annotation           (may be .gz)
        SRR1234567/SRR1234567.sra
        SRR1234568/SRR1234568.sra
    ...
```

- **One `.txt` file per species** in `SRR_to_run/`. The file name (without `.txt`) must match the species folder name in `transcriptomes/`.
- **Each line** of the `.txt` file holds a run accession followed by a free-text condition label, which becomes the column name in the TPM table:
  ```
  SRR1234567 adult female rep1
  SRR1234568 adult male rep1
  ```
- **Exactly one genome and one annotation per species folder.** These can be plain files, `.gz` files (e.g. from WormBase ParaSite), or an NCBI Datasets `.zip` archive.
- **`.sra` files** can be downloaded with the SRA Toolkit, which creates the expected `<SRR_ID>/<SRR_ID>.sra` layout:
  ```bash
  prefetch SRR1234567 -O transcriptomes/Necator_americanus/
  ```
- **Short-read data must be paired-end.** Species with long-read data are listed in `LONG_READ_SPECIES` at the top of the script (default: `Ancylostoma_duodenale`).

## Usage

```bash
python rnaseq_tpm_pipeline.py [METADATA_DIR] [THREADS]
```

For example, to run with the default metadata folder and 32 threads:

```bash
python rnaseq_tpm_pipeline.py SRR_to_run 32
```

Both arguments are optional (defaults: `SRR_to_run`, `32`). `THREADS` applies to all tools except Trimmomatic, which uses 4. Further settings (memory limits, trimming parameters) are in the `PipelineConfiguration` class at the top of the script.

The pipeline can be re-run safely: runs whose BAM file already exists, and species whose results already exist, are skipped.

## Output

| Location | Content |
|---|---|
| `transcriptomes/<Species>/<Species>_03_alignment/` | Sorted and indexed BAM per run, merged BAM, HISAT2/minimap2 alignment summaries |
| `transcriptomes/<Species>/<Species>_05_counts/` | featureCounts table, `<Species>_TPM_Results.csv` / `.xlsx` |
| `All_Species_TPM_Summary_<timestamp>.xlsx` | All TPM tables, one sheet per species |
| `analysis_RNAseq_pipeline_<timestamp>.log` | Full run log, including every command and all software versions |
| `RNAseq_pipeline_errors.csv` | One line per failed command (only created if something fails) |

Intermediate FASTQ files, the HISAT2 index and FastQC reports are deleted once a species has finished, to save disk space.


## License

This project is licensed under the MIT License.
