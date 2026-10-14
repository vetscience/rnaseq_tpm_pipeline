"""
RNA-seq quantification pipeline: SRA runs -> trimmed reads -> genome alignment
-> gene-level counts -> TPM matrices (one per species + one combined Excel workbook).

Workflow per species
--------------------
 1. Reference preparation : unpack genome/annotation (.zip/.gz); GFF3 -> GTF (gffread)
 2. Genome index          : HISAT2 (hisat2-build); not built for the long-read species
 3. Read extraction       : fasterq-dump (SRA Toolkit)
 4. Raw-read QC           : FastQC
 5. Trimming              : Trimmomatic PE (TruSeq3-PE adapters); short-read data must be paired-end
 6. Alignment             : HISAT2 (all short-read species) or minimap2 -ax splice:hq
                            (long reads, Ancylostoma duodenale only); samtools sort + index
 7. Trimmed-read QC       : FastQC
 8. Merging               : samtools merge of all runs of a species ("Merged_Total")
 9. Quantification        : featureCounts (exon features, gene_id, unstranded)
10. Normalisation         : TPM = RPK_i / sum_j(RPK_j) * 1e6, with RPK = count / (length / 1000)
                            and length = featureCounts "Length" column (non-overlapping exon bp)

Expected input layout (relative to this script)
-----------------------------------------------
    SRR_to_run/<Species_name>.txt        one line per run: "<SRR_ID> <condition label ...>"
    transcriptomes/<Species_name>/       exactly one genome and one annotation, plain or .gz,
                                         or an NCBI Datasets .zip:
                                            genome:     *genomic.fna | *genomic.fa | *_genome
                                            annotation: *.gtf | *.gff3 | *.gff | *_annotation (each may be .gz)
    transcriptomes/<Species_name>/<SRR_ID>/<SRR_ID>.sra
    TruSeq3-PE.fa                        Trimmomatic adapter file

Outputs
-------
    transcriptomes/<Species>/<Species>_03_alignment/  sorted + indexed BAMs per run, merged BAM,
                                                      HISAT2/minimap2 alignment summaries
    transcriptomes/<Species>/<Species>_05_counts/     featureCounts table, <Species>_TPM_Results.csv/.xlsx
    All_Species_TPM_Summary_<timestamp>.xlsx          one sheet per species
    analysis_RNAseq_pipeline_<timestamp>.log          full run log, including software versions
    RNAseq_pipeline_errors.csv                        one line per failed command (if any)

Usage
-----
        python rnaseq_tpm_pipeline.py [METADATA_DIR] [THREADS]
        e.g. python rnaseq_tpm_pipeline.py SRR_to_run 32

Requirements: Python >= 3.9, pandas, openpyxl; SRA Toolkit, FastQC, Trimmomatic, HISAT2,
minimap2, samtools, Subread (featureCounts), gffread.

Citation: Kaethner, M., Young, N. D., Sleebs, B. E., Lee, E. F., Fairlie, W. D., & Gasser, R. B. (2026).
Comparative reconstruction of intrinsic apoptosis in clade V nematodes reveals a conserved execution
module and heterogeneous regulatory architecture. Computational and Structural Biotechnology Journal.
https://doi.org/10.34133/csbj.0239
"""

from __future__ import annotations

import argparse
import gzip
import os
import platform
import shlex
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import pandas as pd

TIMEZONE = ZoneInfo("Australia/Melbourne")
def timestamp() -> str:
    """Current time in Melbourne as YYYYmmdd_HHMM (used for log lines and output names)."""
    return datetime.now(TIMEZONE).strftime("%Y%m%d_%H%M")

class PipelineConfiguration:
    """Central configuration of tools, resources and file names."""

    # --- Tools and reference data ---
    TRIMMOMATIC_ADAPTER_FILE = "TruSeq3-PE.fa"
    SRA_DUMP_EXE = "fasterq-dump"

    # --- Folders ---
    INPUT_METADATA_DIR = "SRR_to_run"

    # --- Performance ---
    THREADS = 32                    # all tools except Trimmomatic
    TRIMMOMATIC_THREADS = 4
    FASTERQ_MEM_LIMIT = "150G"
    SAMTOOLS_MEM_PER_THREAD = "1G"

    # --- Long-read species (minimap2, featureCounts without -p). All other species: HISAT2. ---
    LONG_READ_SPECIES: list[str] = ["Ancylostoma_duodenale"]

    # --- Trimmomatic settings (identical to the original analysis) ---
    TRIM_ILLUMINACLIP = "2:30:10"
    TRIM_STEPS = ["LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15"]
    TRIM_MINLEN = 36

    # --- Outputs ---
    _now = timestamp()
    FINAL_EXCEL = f"All_Species_TPM_Summary_{_now}.xlsx"
    LOG_FILE = f"analysis_RNAseq_pipeline_{_now}.log"
    ERROR_FILE = "RNAseq_pipeline_errors.csv"

PC = PipelineConfiguration()


def log(msg):
    """Print a timestamped message and append it to the run log."""
    line = f"[{timestamp()}] {msg}"
    print(line, flush=True)
    with open(PC.LOG_FILE, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")

def record_error(species_name: str, sample: str, step: str, detail: str) -> None:
    """Append one line to the error table (species, sample, step, detail)."""
    detail = detail.replace(",", ";").replace("\n", " ")
    with open(PC.ERROR_FILE, "a", encoding="utf-8") as fh:
        fh.write(f"{species_name},{sample},{step},{detail}\n")

def run_cmd(cmd: list[str], species_name: str, sample: str = "-") -> None:
    """Run one external command; stdout/stderr go to the run log. Raise RuntimeError on failure."""
    cmd = [str(c) for c in cmd]
    log(f"Executing: {shlex.join(cmd)} for {species_name}")
    with open(PC.LOG_FILE, "a", encoding="utf-8") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=fh)
    if result.returncode != 0:
        msg = f"[ERROR] Command failed with exit code {result.returncode}: {shlex.join(cmd)}"
        log(msg)
        record_error(species_name, sample, cmd[0], f"exit={result.returncode} {shlex.join(cmd)}")
        raise RuntimeError(msg)
def run_piped_commands(commands: list[list[str]], species_name: str, sample: str,
                       aligner_log: Path) -> None:
    """Run commands connected by pipes, e.g. `hisat2 ... | samtools view ... | samtools sort ...`.

    commands[0] is the aligner; its output is passed to commands[1], whose output is passed to
    commands[2], and so on. The aligner's messages (for HISAT2: the alignment-rate summary) are
    saved to `aligner_log`. `set -o pipefail` makes the whole pipe fail if ANY of the commands
    fails; without it only the exit code of the last one (samtools sort) would be checked.
    """
    shell_parts = [shlex.join([str(arg) for arg in cmd]) for cmd in commands]
    shell_parts[0] += f" 2> {shlex.quote(str(aligner_log))}"
    shell_cmd = "set -o pipefail; " + " | ".join(shell_parts)
    log(f"Executing: {shell_cmd} for {species_name}")
    result = subprocess.run(shell_cmd, shell=True, executable="/bin/bash")
    if result.returncode != 0:
        msg = f"[ERROR] Command failed with exit code {result.returncode}: {shell_cmd}"
        log(msg)
        record_error(species_name, sample, commands[0][0], f"exit={result.returncode} {shell_cmd}")
        raise RuntimeError(msg)

SOFTWARE_VERSION_COMMANDS = {
    "fasterq-dump": [PC.SRA_DUMP_EXE, "--version"],
    "FastQC": ["fastqc", "--version"],
    "Trimmomatic": ["trimmomatic", "-version"],
    "HISAT2": ["hisat2", "--version"],
    "minimap2": ["minimap2", "--version"],
    "samtools": ["samtools", "--version"],
    "featureCounts": ["featureCounts", "-v"],
    "gffread": ["gffread", "--version"],
}


def log_software_versions() -> None:
    """Write Python, package and external tool versions to the log."""
    log(f"Python {platform.python_version()} on {platform.platform()}")
    log(f"pandas {pd.__version__}")
    try:
        import openpyxl
        log(f"openpyxl {openpyxl.__version__}")
    except ImportError:
        log("openpyxl: NOT INSTALLED (needed for Excel output)")
    for name, cmd in SOFTWARE_VERSION_COMMANDS.items():
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, errors="replace", timeout=120)
            lines = [ln.strip() for ln in (res.stdout + res.stderr).splitlines() if ln.strip()]
            log(f"{name}: {lines[0] if lines else 'version not reported'}")
        except FileNotFoundError:
            log(f"{name}: NOT FOUND in PATH")
        except (subprocess.TimeoutExpired, OSError) as e:
            log(f"{name}: version check failed ({e})")

def decompress_reference_files(species_dir: Path, species_name: str) -> tuple[Path | None, Path | None]:
    """Locate (and if needed unzip/gunzip) the genome FASTA and annotation of one species."""
    log(f"-- Checking genome and annotation files for {species_name}...")

    search_patterns = {
        "genome": ["*genomic.fna", "*genomic.fa", "*_genome"],
        "annotation": ["*.gtf", "*.gff3", "*.gff", "*_annotation"],
    }

    def scan() -> dict[str, Path | None]:
        """Find uncompressed genome/annotation files (one of each per species folder)."""
        res: dict[str, Path | None] = {"genome": None, "annotation": None}
        for key, patterns in search_patterns.items():
            for pattern in patterns:
                files = list(species_dir.glob(pattern))
                if files:
                    res[key] = files[0]
                    break
        return res

    found_paths = scan()

    # NCBI Datasets (.zip archives)
    if not found_paths["genome"] or not found_paths["annotation"]:
        for zf in sorted(species_dir.glob("*.zip")):
            log(f"Extracting zip: {zf.name}...")
            with zipfile.ZipFile(zf, 'r') as zip_ref:
                zip_ref.extractall(species_dir)
            ncbi_path = species_dir / "ncbi_dataset"
            if ncbi_path.exists():
                for f in ncbi_path.rglob("*"):
                    if f.is_file():
                        dest = species_dir / f.name
                        if not dest.exists():
                            shutil.move(str(f), str(dest))
                shutil.rmtree(ncbi_path)
        found_paths = scan()

    # WormBase ParaSite Datasets (.gz files)
    for key, patterns in search_patterns.items():
        if not found_paths[key]:
            for pattern in patterns:
                gz_files = list(species_dir.glob(pattern + ".gz"))
                if gz_files:
                    gz_path = gz_files[0]
                    uncompressed_path = gz_path.with_suffix('')
                    if not uncompressed_path.exists():
                        log(f"Decompressing {gz_path.name}...")
                        with gzip.open(gz_path, 'rb') as f_in, open(uncompressed_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)

                    found_paths[key] = uncompressed_path
                    break

    return found_paths["genome"], found_paths["annotation"]


def convert_gff3_to_gtf(species_dir: Path, species_name: str) -> Path | None:
    """Return a GTF for the species, converting from GFF3 with gffread if needed."""
    log("-- Checking whether a GTF annotation is available --")

    gtf_files = list(species_dir.glob("*.gtf"))
    if gtf_files:
        gtf_path = gtf_files[0]
        if gtf_path.stat().st_size > 0:
            log(f"Valid GTF found ({gtf_path.name}, {gtf_path.stat().st_size} bytes). Skipping.")
            return gtf_path
        log(f"Found empty GTF ({gtf_path.name}). Re-converting.")

    gff3_files = [f for f in species_dir.glob("*")
                  if f.name.endswith((".gff3", "annotation", ".gff")) and not f.name.endswith(".gtf")]
    if not gff3_files:
        log(f"[WARNING] No GTF or GFF3 annotation found in {species_dir}")
        return None
    gff3_path = gff3_files[0]
    log(f"Found GFF3 file to convert: {gff3_path.name}")

    # --keep-genes and -F make WormBase ParaSite GFF3 convert to the same GTF layout as NCBI GTFs.
    gtf_path = gff3_path.with_suffix(".gtf")
    log(f"Converting {gff3_path.name} to {gtf_path.name} using gffread...")
    try:
        run_cmd(["gffread", "-T", str(gff3_path), "--keep-genes", "-F", "-o", str(gtf_path)], species_name)
        log(f"GFF3 -> GTF conversion for {species_name} complete: {gtf_path.name}")
        return gtf_path
    except RuntimeError as e:
        log(f"[ERROR] gffread conversion for {species_name} failed: {e}")
        return None

def get_metadata(txt_path: Path) -> dict[str, str]:
    """Read '<SRR_ID> <condition ...>' lines into {SRR_ID: 'condition (SRR_ID)'}."""
    mapping = {}
    with open(txt_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                mapping[parts[0]] = f"{' '.join(parts[1:])} ({parts[0]})"
            elif len(parts) == 1:
                mapping[parts[0]] = parts[0]
    return mapping

def build_hisat2_index(ref_genome: Path, index_prefix: Path, species_name: str) -> None:
    """Build the HISAT2 genome index unless it already exists."""
    if all(Path(f"{index_prefix}.{i}.ht2").exists() for i in range(1, 9)):
        log(f"HISAT2 index for {species_name} already exists - skipping this step.")
        return
    log(f"-- STEP 2: BUILD HISAT2 INDEX for {species_name}")
    index_prefix.parent.mkdir(parents=True, exist_ok=True)
    run_cmd(["hisat2-build", "-p", str(PC.THREADS), str(ref_genome), str(index_prefix)], species_name)
    log("HISAT2 index built successfully.")

def run_fasterq_dump(srr_id: str, sra_file: Path, output_dir: Path, species_name: str,
                     temp_dir: Path) -> tuple[Path, Path | None]:
    """Convert one .sra run to FASTQ. Returns (R1, R2) for paired data, (FASTQ, None) for long reads."""
    log(f"-- STEP 3: FASTERQ-DUMP (SRA -> FASTQ) for {srr_id}")

    single_fastq = output_dir / f"{srr_id}.fastq"
    r1 = output_dir / f"{srr_id}_1.fastq"
    r2 = output_dir / f"{srr_id}_2.fastq"

    if single_fastq.exists():
        log(f"Found existing single-end FASTQ for {srr_id}.")
        return single_fastq, None

    # Incomplete paired output from an interrupted run
    if (r1.exists() and not r2.exists()) or (r2.exists() and not r1.exists()):
        log(f"[FIX] Incomplete paired FASTQ detected for {srr_id}, cleaning up.")
        if r1.exists(): r1.unlink()
        if r2.exists(): r2.unlink()

    # Empty files (silent fasterq-dump failure)
    for f in [r1, r2]:
        if f.exists() and f.stat().st_size == 0:
            log(f"[FIX] Empty FASTQ detected: {f.name}, deleting.")
            f.unlink()

    if r1.exists() and r2.exists():
        log(f"FASTQ files for {srr_id} ({species_name}) already exist. Skipping extraction.")
        return r1, r2

    if not sra_file.exists():
        raise FileNotFoundError(f"SRA file not found at {sra_file}")

    run_cmd([PC.SRA_DUMP_EXE, str(sra_file),
             "--split-files",
             "--threads", str(PC.THREADS),
             "--mem", PC.FASTERQ_MEM_LIMIT,
             "--outdir", str(output_dir),
             "--temp", str(temp_dir)], species_name, srr_id)

    if r1.exists() and r2.exists():
        return r1, r2
    elif r1.exists():
        return r1, None
    elif single_fastq.exists():
        return single_fastq, None
    else:
        raise RuntimeError(f"fasterq-dump finished but no FASTQ files found for {srr_id}.")


def run_fastqc(fastq_files: list[Path | None], out_dir: Path, species_name: str, label: str) -> None:
    """Run FastQC on each existing FASTQ unless its report already exists."""
    log(f"-- FASTQC ({label} reads) for {species_name}")
    for fastq in (Path(f) for f in fastq_files if f is not None and Path(f).exists()):
        if (out_dir / f"{fastq.stem}_fastqc.html").exists():
            log(f"FastQC report for {fastq.name} already exists - skipping.")
            continue
        run_cmd(["fastqc", "-t", str(PC.THREADS), "-o", str(out_dir), str(fastq)], species_name)



def run_trimmomatic(srr_id: str, fastq_r1: Path, fastq_r2: Path, trim_dir: Path,
                    species_name: str) -> tuple[Path, Path]:
    """Adapter- and quality-trim one paired-end run. Returns (trimmed R1, trimmed R2)."""
    log(f"-- STEP 5: TRIMMOMATIC for {srr_id} ({species_name})")

    trimmed_r1_paired = trim_dir / f'{srr_id}_R1_paired.fastq'
    trimmed_r1_unpaired =trim_dir / f'{srr_id}_R1_unpaired.fastq'
    trimmed_r2_paired = trim_dir / f'{srr_id}_R2_paired.fastq'
    trimmed_r2_unpaired = trim_dir / f'{srr_id}_R2_unpaired.fastq'

    if trimmed_r1_paired.exists() and trimmed_r2_paired.exists():
        log(f"Trimming skipped for {srr_id}: output files already exist.")
        return trimmed_r1_paired, trimmed_r2_paired

    cmd = (["trimmomatic", "PE", "-phred33", "-threads", str(PC.TRIMMOMATIC_THREADS),
            str(fastq_r1), str(fastq_r2),
            str(trimmed_r1_paired), str(trimmed_r1_unpaired),
            str(trimmed_r2_paired), str(trimmed_r2_unpaired),
            f"ILLUMINACLIP:{PC.TRIMMOMATIC_ADAPTER_FILE}:{PC.TRIM_ILLUMINACLIP}"]
           + PC.TRIM_STEPS + [f"MINLEN:{PC.TRIM_MINLEN}"])

    run_cmd(cmd, species_name, srr_id)
    log(f"Trimmomatic completed successfully for {srr_id}")
    return trimmed_r1_paired, trimmed_r2_paired

def _remove_partial_bam(out_bam: Path, temp_prefix: str | None = None) -> None:
    """Delete a BAM (and its index/temp files) left behind by a failed alignment."""
    for f in [out_bam, Path(str(out_bam) + ".bai")]:
        f.unlink(missing_ok=True)
    if temp_prefix:
        for f in Path(temp_prefix).parent.glob(Path(temp_prefix).name + "*"):
            f.unlink(missing_ok=True)

def run_alignment(trimmed_r1: Path, trimmed_r2: Path, srr_id: str, align_dir: Path,
                  index_prefix: Path, species_name: str) -> Path:
    """Align trimmed paired-end reads with HISAT2; sort and index the BAM with samtools."""
    out_bam_sorted = align_dir / f"{srr_id}.sorted.bam"
    out_bai_index = align_dir / f"{srr_id}.sorted.bam.bai"
    temp_prefix = str(out_bam_sorted) + ".tmp"
    log(f"-- STEP 6: HISAT2 ALIGNMENT + SAMTOOLS SORT for {srr_id} ({species_name})")

    if out_bam_sorted.exists() and out_bai_index.exists():
        log(f"Skipping HISAT2 for {srr_id}: {out_bam_sorted.name} and index already exist.")
        return out_bam_sorted

    hisat2_cmd = ["hisat2", "-x", str(index_prefix),
                  "-1", str(trimmed_r1), "-2", str(trimmed_r2),
                  "--threads", str(PC.THREADS)]
    view_cmd = ["samtools", "view", "-@", str(PC.THREADS), "-b"]
    sort_cmd = ["samtools", "sort", "-@", str(PC.THREADS), "-T", temp_prefix,
                "-o", str(out_bam_sorted), "-m", PC.SAMTOOLS_MEM_PER_THREAD]

    try:
        run_piped_commands([hisat2_cmd, view_cmd, sort_cmd], species_name, srr_id,
                           aligner_log=align_dir / f"{srr_id}.hisat2_summary.txt")
    except RuntimeError:
        _remove_partial_bam(out_bam_sorted, temp_prefix)
        raise

    run_cmd(["samtools", "index", str(out_bam_sorted)], species_name, srr_id)
    log(f"Alignment + sorting complete: {out_bam_sorted}")
    return out_bam_sorted


def run_samtools_merge(align_dir: Path, sorted_bam_files: list[str], species_name: str) -> Path:
    """Merge all per-run BAMs of a species into one BAM (reported as 'Merged_Total')."""
    log(f"-- STEP 8: SAMTOOLS MERGE for {species_name}")
    merged_bam = align_dir / f'{species_name}_merged.bam'

    if merged_bam.exists():
        log(f"Merged BAM already exists: {merged_bam.name}. Skipping.")
        return merged_bam

    run_cmd(["samtools", "merge", "-@", str(PC.THREADS), str(merged_bam)] + sorted_bam_files, species_name)
    run_cmd(["samtools", "index", str(merged_bam)], species_name)
    return merged_bam

def run_minimap2_longread(fastq_file: Path, ref_genome: Path, align_dir: Path, species_name: str) -> Path:
    """Spliced long-read alignment with minimap2 (-ax splice:hq -uf); sort and index the BAM."""
    log(f"-- STEP 6: MINIMAP2 LONG-READ ALIGNMENT for {species_name}")
    run_name = Path(fastq_file).stem
    out_bam_sorted = align_dir / f"{run_name}.sorted.bam"
    out_bai = align_dir / f"{run_name}.sorted.bam.bai"

    if out_bam_sorted.exists() and out_bai.exists():
        log(f"Long-read BAM already exists: {out_bam_sorted.name}")
        return out_bam_sorted

    minimap_cmd = ["minimap2", "-t", str(PC.THREADS), "-ax", "splice:hq", "-uf",
                   str(ref_genome), str(fastq_file)]
    view_cmd = ["samtools", "view", "-@", str(PC.THREADS), "-b"]
    sort_cmd = ["samtools", "sort", "-@", str(PC.THREADS), "-o", str(out_bam_sorted),
                "-m", PC.SAMTOOLS_MEM_PER_THREAD]
    try:
        run_piped_commands([minimap_cmd, view_cmd, sort_cmd], species_name, run_name,
                           aligner_log=align_dir / f"{run_name}.minimap2_log.txt")
    except RuntimeError:
        _remove_partial_bam(out_bam_sorted)
        raise

    run_cmd(["samtools", "index", str(out_bam_sorted)], species_name, run_name)
    return out_bam_sorted

def run_featurecounts(sorted_bam_files: list[str], annotation_gtf: Path, count_dir: Path,
                      species_name: str) -> Path:
    """Count reads per gene (exon features summarised by gene_id)."""
    log(f"-- STEP 9: FEATURECOUNTS for {species_name}")
    out_count_file = count_dir / f"gene_counts_{species_name}.txt"
    if out_count_file.exists():
        log(f"Skipping featureCounts: {out_count_file} already exists.")
        return out_count_file

    # Paired-end counting (-p -B -C) for all short-read species; plain read counting for long reads.
    is_long_read = species_name in PC.LONG_READ_SPECIES
    paired_flags = [] if is_long_read else ["-p", "-B", "-C"]
    if is_long_read:
        log("Long-read data: featureCounts run without -p.")

    cmd = (["featureCounts"] + paired_flags
           + ["-t", "exon", "-g", "gene_id", "-s", "0",
              "-a", str(annotation_gtf), "-o", str(out_count_file),
              "-T", str(PC.THREADS), "-F", "GTF"]
           + sorted_bam_files)
    run_cmd(cmd, species_name)
    return out_count_file
FEATURECOUNTS_ANNOTATION_COLUMNS = ["Chr", "Start", "End", "Strand", "Length"]

def calculate_tpm_and_export_all(count_file: Path, srr_mapping: dict[str, str], count_dir: Path,
                                 species_name: str) -> pd.DataFrame:
    """Compute TPM from a featureCounts table and write CSV + XLSX.

    RPK_i = count_i / (length_i / 1000);  TPM_i = RPK_i / sum_j(RPK_j) * 1e6,
    with length = featureCounts 'Length' (number of non-overlapping exonic bases of the gene).
    """
    log(f"-- STEP 10: TPM CALCULATION for {species_name}")

    # First line of the featureCounts output is the program/command comment.
    df = pd.read_csv(count_file, sep="\t", skiprows=1, index_col=0)
    counts_df = df.drop(columns=FEATURECOUNTS_ANNOTATION_COLUMNS)
    gene_lengths = df["Length"]

    rpk = counts_df.div(gene_lengths.values / 1000, axis=0)   # reads per kilobase (RPK)
    # Divide each sample's RPK by (its RPK total / 1e6), i.e. TPM = RPK / total * 1e6,
    # so every sample's TPM values add up to 1,000,000.
    tpm_df = rpk.div(rpk.sum(axis=0) / 1e6, axis=1)

    # Column names: BAM path -> SRR ID -> "condition (SRR ID)" from the metadata file
    tpm_df.columns = [Path(c).name.replace(".sorted.bam", "").replace(".bam", "") for c in tpm_df.columns]
    column_labels = dict(srr_mapping)
    column_labels[f"{species_name}_merged"] = "Merged_Total"
    tpm_df = tpm_df.rename(columns=column_labels)

    out_base = count_dir / f"{species_name}_TPM_Results"
    tpm_df.to_csv(f"{out_base}.csv")
    tpm_df.to_excel(f"{out_base}.xlsx")
    log(f"TPM calculation complete. Results saved to {out_base}.csv/.xlsx")
    return tpm_df

def process_species(species_name: str, srr_mapping: dict[str, str]) -> pd.DataFrame | None:
    """Run the complete workflow for one species and return its TPM table (None on failure)."""
    log(f"\n{'#' * 60}\nSTARTING ANALYSIS: {species_name}\n{'#' * 60}")

    base_dir = Path.cwd() / "transcriptomes" / species_name
    qc_init_dir = base_dir / f"{species_name}_01_initial_qc"
    trim_dir = base_dir / f"{species_name}_02_trimmed_fastq"
    align_dir = base_dir / f"{species_name}_03_alignment"
    qc_trim_dir = base_dir / f"{species_name}_04_trimmed_qc"
    count_dir = base_dir / f"{species_name}_05_counts"
    temp_dir = base_dir / f"{species_name}_06_temp"
    for d in (qc_init_dir, trim_dir, align_dir, qc_trim_dir, count_dir, temp_dir):
        d.mkdir(parents=True, exist_ok=True)


    tpm_csv = count_dir / f"{species_name}_TPM_Results.csv"
    if tpm_csv.exists() and (align_dir / f"{species_name}_merged.bam").exists():
        log(f"SKIPPING {species_name}: final results already exist.")
        return pd.read_csv(tpm_csv, index_col=0)

    is_long_read = species_name in PC.LONG_READ_SPECIES
    index_prefix = base_dir / "hisat2_index" / f"{species_name}_idx"

    # --- STEP 1-2: reference preparation and genome index ---
    try:
        genome, _ = decompress_reference_files(base_dir, species_name)
        annotation = convert_gff3_to_gtf(base_dir, species_name)
        if not genome or not annotation:
            log(f"SKIPPING: missing genome or annotation in {base_dir}.")
            return None
        if not is_long_read:
            build_hisat2_index(genome, index_prefix, species_name)
    except Exception as e:
        log(f"[ERROR] Could not prepare reference files for {species_name}: {e}")
        return None

    sra_files = sorted(base_dir.glob("**/*.sra"))
    if not sra_files:
        log(f"[ERROR] No .sra files found for {species_name} in {base_dir}.")
        return None

    all_sorted_bams: list[str] = []
    failed_samples: list[str] = []

    # --- STEP 3-7: per-run processing ---
    for sra_path in sra_files:
        srr_id = sra_path.stem
        fastq_r1 = fastq_r2 = trimmed_r1 = trimmed_r2 = None
        out_bam = align_dir / f'{srr_id}.sorted.bam'
        if out_bam.exists():
            log(f"Found existing BAM for {srr_id}, adding to list.")
            all_sorted_bams.append(str(out_bam))
            continue
        log(f"-- PROCESSING RUN {srr_id} for {species_name}")
        try:
            fastq_r1, fastq_r2 = run_fasterq_dump(srr_id, sra_path, base_dir, species_name, temp_dir)

            if is_long_read:
                log(f"{species_name} is a long-read species: aligning with minimap2.")
                out_bam = run_minimap2_longread(fastq_r1, genome, align_dir, species_name)
            else:
                if fastq_r2 is None:
                    raise RuntimeError(f"{srr_id} is not paired-end; this pipeline expects "
                                       f"paired-end short reads for {species_name}.")
                run_fastqc([fastq_r1, fastq_r2], qc_init_dir, species_name, "raw")
                trimmed_r1, trimmed_r2 = run_trimmomatic(srr_id, fastq_r1, fastq_r2, trim_dir, species_name)
                out_bam = run_alignment(trimmed_r1, trimmed_r2, srr_id, align_dir, index_prefix, species_name)
                run_fastqc([trimmed_r1, trimmed_r2], qc_trim_dir, species_name, "trimmed")

            all_sorted_bams.append(str(out_bam))

            # Remove raw and trimmed FASTQ files once the BAM is safely written
            if out_bam.exists() and out_bam.stat().st_size > 0:
                log(f"Success for {srr_id}. Removing intermediate FASTQ files...")
                for f in (fastq_r1, fastq_r2, trimmed_r1, trimmed_r2):
                    if f is not None and Path(f).exists():
                        Path(f).unlink()
                        log(f"Deleted intermediate: {Path(f).name}")



        except (RuntimeError, OSError) as e:
            failed_samples.append(srr_id)
            record_error(species_name, srr_id, "sample", str(e))
            log(f"[ERROR] Sample {srr_id} ({species_name}) failed: {e}. Continuing with next sample.")

    if failed_samples:
        log(f"[WARNING] {len(failed_samples)} run(s) failed for {species_name} and are NOT in the "
            f"TPM table: {', '.join(failed_samples)}")

    if not all_sorted_bams:
        return None

    # --- STEP 8-10: merge, count, normalise ---
    merged_bam = run_samtools_merge(align_dir, all_sorted_bams, species_name)
    count_file = run_featurecounts(all_sorted_bams + [str(merged_bam)], annotation, count_dir, species_name)
    tpm_df = calculate_tpm_and_export_all(count_file, srr_mapping, count_dir, species_name)

    # Keep SRA runs, alignments and counts; remove everything that can be regenerated
    log(f"Removing temporary folders for {species_name}...")
    for d in (base_dir / "hisat2_index", trim_dir, qc_init_dir, qc_trim_dir, temp_dir):
        shutil.rmtree(d, ignore_errors=True)

    log(f"\n{'#' * 60}\nCOMPLETED RNA-SEQ PIPELINE FOR: {species_name}\n{'#' * 60}")
    return tpm_df

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Read command-line options (defaults are the values used in the paper)."""
    parser = argparse.ArgumentParser(
        description="RNA-seq pipeline: SRA -> FastQC/Trimmomatic -> HISAT2/minimap2 -> featureCounts -> TPM.")
    parser.add_argument("metadata_dir", nargs="?", default=PC.INPUT_METADATA_DIR,
                        help="Folder with one <Species_name>.txt file per species (default: %(default)s)")
    parser.add_argument("threads", nargs="?", type=int, default=PC.THREADS,
                        help="Threads for all tools except Trimmomatic (default: %(default)s)")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the pipeline for every species listed in the metadata folder."""
    # All relative paths (metadata, transcriptomes/, adapter file) are resolved from the script folder.
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    args = parse_args(argv)
    PC.THREADS = args.threads

    input_path = Path(args.metadata_dir)
    if not input_path.exists():
        print(f"[ERROR] Metadata directory not found: {input_path}", file=sys.stderr)
        return 1

    metadata_files = sorted(input_path.glob("*.txt"))
    if not metadata_files:
        print(f"[WARNING] No .txt files found in {input_path}", file=sys.stderr)
        return 1

    log(f"Pipeline started. Command: {shlex.join(sys.argv)}")
    log_software_versions()

    all_tpm_results: dict[str, pd.DataFrame] = {}
    for txt_file in metadata_files:
        species_name = txt_file.stem
        mapping = get_metadata(txt_file)
        try:
            tpm_data = process_species(species_name, mapping)
            if tpm_data is not None:
                all_tpm_results[species_name] = tpm_data
        except Exception as e:
            log(f"CRITICAL FAILURE for species {species_name}: {e}")
            record_error(species_name, "-", "species", str(e))

    log(f"Results produced for {len(all_tpm_results)} of {len(metadata_files)} species.")

    if not all_tpm_results:
        log("[WARNING] No TPM data was generated. Combined Excel file skipped.")
        return 1

    log("=" * 50)
    log(f"Writing combined Excel workbook: {PC.FINAL_EXCEL}")
    try:
        with pd.ExcelWriter(PC.FINAL_EXCEL, engine="openpyxl") as writer:
            for species, df in all_tpm_results.items():
                df.to_excel(writer, sheet_name=species[:31])  # Excel sheet-name limit: 31 characters
        log("Combined Excel workbook saved successfully.")
    except Exception as e:
        log(f"[ERROR] Failed to save combined Excel workbook: {e}")
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())