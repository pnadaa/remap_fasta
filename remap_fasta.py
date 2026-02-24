#!/usr/bin/env python3
"""
remap_fasta_coords.py

For each sequence in a multifasta (headers: accessionID_start-end), aligns the
sequence to its reference genome region via blastn and outputs an updated
multifasta with precise, remapped start/end coordinates.

Dependencies: BLAST+ (blastn, blastdbcmd), Biopython
"""

import argparse
import os
import subprocess
import sys
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------
# Header parsing
# ---------------------------------------------------------------------------

def parse_header(header: str) -> tuple[str, int, int]:
    """
    Parse a header of the form  accessionID_start-end.
    Uses rsplit('_', 1) so accession IDs containing underscores are handled
    correctly.  Returns (accession, start, end) as 1-based coordinates.
    """
    try:
        accession, coord_str = header.rsplit("_", 1)
        start_str, end_str = coord_str.split("-", 1)
        return accession, int(start_str), int(end_str)
    except (ValueError, AttributeError) as exc:
        raise ValueError(
            f"Cannot parse header '{header}': expected accessionID_start-end"
        ) from exc


# ---------------------------------------------------------------------------
# Reference region extraction
# ---------------------------------------------------------------------------

def extract_region(db: str, accession: str, start: int, end: int) -> str:
    """
    Use blastdbcmd to pull the genomic region [start, end] (1-based, inclusive)
    from the local BLAST database.  Returns a FASTA string.
    """
    cmd = [
        "blastdbcmd",
        "-db",    db,
        "-entry", accession,
        "-range", f"{start}-{end}",
        "-outfmt", "%f",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not result.stdout.strip():
        raise RuntimeError(
            f"blastdbcmd failed for {accession}:{start}-{end}\n"
            f"  STDERR: {result.stderr.strip()}"
        )
    return result.stdout


# ---------------------------------------------------------------------------
# BLAST alignment
# ---------------------------------------------------------------------------

def run_blastn(
    query_fasta: str,
    subject_fasta: str,
    task: str = "blastn",
    perc_identity: float = 80.0,
    max_hsps: int = 1,
) -> str:
    """
    Align query_fasta (FASTA string) against subject_fasta (FASTA string)
    using blastn.  Returns the raw tabular (outfmt 6) output string.
    """
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as qf:
        qf.write(query_fasta)
        query_file = qf.name

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as sf:
        sf.write(subject_fasta)
        subject_file = sf.name

    try:
        cmd = [
            "blastn",
            "-query",        query_file,
            "-subject",      subject_file,
            "-task",         task,
            "-perc_identity", str(perc_identity),
            "-max_hsps",     str(max_hsps),
            "-outfmt",       "6 qseqid sseqid pident length mismatch gapopen "
                             "qstart qend sstart send evalue bitscore",
            "-dust",         "no",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"blastn failed:\n  {result.stderr.strip()}")
        return result.stdout
    finally:
        os.unlink(query_file)
        os.unlink(subject_file)


# ---------------------------------------------------------------------------
# BLAST result parsing
# ---------------------------------------------------------------------------

def best_hit(blast_output: str) -> tuple[int, int] | None:
    """
    Parse blastn outfmt-6 output and return (sstart, send) for the HSP
    with the highest bitscore.  sstart > send indicates a minus-strand hit.
    Returns None if there are no hits.
    """
    best_coords = None
    best_score  = -1.0

    for line in blast_output.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 12:
            continue
        try:
            sstart   = int(cols[8])
            send     = int(cols[9])
            bitscore = float(cols[11])
        except ValueError:
            continue
        if bitscore > best_score:
            best_score  = bitscore
            best_coords = (sstart, send)

    return best_coords


# ---------------------------------------------------------------------------
# Coordinate conversion
# ---------------------------------------------------------------------------

def relative_to_absolute(
    region_start: int, rel_start: int, rel_end: int
) -> tuple[int, int]:
    """
    Convert subject-relative 1-based coordinates from a blastdbcmd-extracted
    region back to absolute genome coordinates.

    region_start : 1-based start of the extracted region in the genome
    rel_start    : BLAST sstart (1-based, within extracted region)
    rel_end      : BLAST send   (1-based, within extracted region)

    Strand orientation is preserved:
      rel_start < rel_end  →  forward strand
      rel_start > rel_end  →  reverse complement (abs_start > abs_end)
    """
    offset    = region_start - 1
    abs_start = rel_start + offset
    abs_end   = rel_end   + offset
    return abs_start, abs_end


# ---------------------------------------------------------------------------
# FASTA output helper
# ---------------------------------------------------------------------------

def write_fasta(records: list[SeqRecord], path: str) -> None:
    """Write records to path as plain FASTA (no trailing description space)."""
    with open(path, "w") as fh:
        for rec in records:
            fh.write(f">{rec.id}\n{str(rec.seq)}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Remap multifasta coordinates by aligning each sequence to its "
            "reference genome region via blastn.\n"
            "Input headers must follow the format:  accessionID_start-end"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input",    required=True,
        help="Input multifasta (headers: accessionID_start-end)",
    )
    parser.add_argument(
        "-o", "--output",   required=True,
        help="Output multifasta with updated coordinates",
    )
    parser.add_argument(
        "-db", "--database", required=True,
        help="Path to local BLAST nucleotide database (blastdbcmd + blastn)",
    )
    parser.add_argument(
        "--task", default="blastn",
        choices=["blastn", "blastn-short", "dc-megablast", "megablast"],
        help="blastn task (default: blastn; use blastn-short for sequences <50 bp)",
    )
    parser.add_argument(
        "--perc_identity", type=float, default=80.0,
        help="Minimum %%identity for BLAST hits (default: 80.0)",
    )
    parser.add_argument(
        "--keep_on_fail", action="store_true",
        help=(
            "Retain records unchanged in output when alignment fails. "
            "Default: omit failed records and report to stderr."
        ),
    )
    args = parser.parse_args()

    records_out: list[SeqRecord] = []
    n_updated = n_failed = 0

    with open(args.input) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            seq_str = str(record.seq)

            # ── 1. Parse header ─────────────────────────────────────────────
            try:
                accession, reg_start, reg_end = parse_header(record.id)
            except ValueError as exc:
                print(f"[WARNING] {exc} — record skipped/kept.", file=sys.stderr)
                n_failed += 1
                if args.keep_on_fail:
                    records_out.append(record)
                continue

            # ── 2. Extract reference region via blastdbcmd ──────────────────
            try:
                region_fasta = extract_region(
                    args.database, accession, reg_start, reg_end
                )
            except RuntimeError as exc:
                print(f"[WARNING] {exc} — record skipped/kept.", file=sys.stderr)
                n_failed += 1
                if args.keep_on_fail:
                    records_out.append(record)
                continue

            # ── 3. Run blastn ───────────────────────────────────────────────
            query_fasta = f">{record.id}\n{seq_str}\n"
            try:
                blast_out = run_blastn(
                    query_fasta, region_fasta,
                    task=args.task,
                    perc_identity=args.perc_identity,
                )
            except RuntimeError as exc:
                print(
                    f"[WARNING] blastn error for {record.id}: {exc}",
                    file=sys.stderr,
                )
                n_failed += 1
                if args.keep_on_fail:
                    records_out.append(record)
                continue

            # ── 4. Parse best HSP ───────────────────────────────────────────
            hit = best_hit(blast_out)
            if hit is None:
                print(
                    f"[WARNING] No BLAST hit for {record.id} — try lowering "
                    f"--perc_identity or switching to --task blastn-short.",
                    file=sys.stderr,
                )
                n_failed += 1
                if args.keep_on_fail:
                    records_out.append(record)
                continue

            # ── 5. Convert to absolute genome coordinates ───────────────────
            rel_s, rel_e = hit
            abs_s, abs_e = relative_to_absolute(reg_start, rel_s, rel_e)

            new_id = f"{accession}_{abs_s}-{abs_e}"
            print(f"[INFO]  {record.id}  →  {new_id}", file=sys.stderr)

            records_out.append(
                SeqRecord(Seq(seq_str), id=new_id, description="", name="")
            )
            n_updated += 1

    # ── 6. Write output ──────────────────────────────────────────────────────
    write_fasta(records_out, args.output)
    print(
        f"\n[DONE]  {n_updated} records updated, {n_failed} failed. "
        f"Output → {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
