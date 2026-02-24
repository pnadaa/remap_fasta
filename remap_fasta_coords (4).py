#!/usr/bin/env python3
"""
remap_fasta_coords.py

For each sequence in a multifasta (headers: accessionID_start-end), aligns the
sequence to its reference genome region via blastn and outputs an updated
multifasta with precise, remapped start/end coordinates.

Dependencies: BLAST+ (blastn, blastdbcmd), Biopython

Valid blastn reward/penalty -> gapopen/gapextend combinations (NCBI Table D1):
  1/-1 : (3,2) (2,2) (1,2) (0,2) (4,1) (3,1) (2,1) (1,1) (0,1)
  1/-2 : (4,4) (2,4) (0,4) (3,3) (1,3) (6,2) (5,2) (4,2) (2,2) (1,2) (0,2) [default: (2,2)]
  1/-3 : (5,2) (2,2) (1,2) (0,2) (2,1) (1,1)
  2/-3 : (5,2) (2,2) (1,2) (0,2) (3,1) (2,1) (1,1)  [blastn default]
  4/-5 : (12,8) (6,5) (5,5) (4,5) (3,5)
  1/-4 : (1,2) (0,2) (2,1) (1,1)
"""

import argparse
import os
import subprocess
import sys
import tempfile
from functools import partial
from multiprocessing import Pool
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
        "-db",     db,
        "-entry",  accession,
        "-range",  f"{start}-{end}",
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
    task: str          = "blastn",
    perc_identity: float = 80.0,
    reward: int        = 1,
    penalty: int       = -2,
    gapopen: int       = 2,
    gapextend: int     = 2,
) -> str:
    """
    Align query_fasta (FASTA string) against subject_fasta (FASTA string)
    using blastn.  Returns the raw tabular (outfmt 6) output string.

    Scoring defaults (reward=1, penalty=-2, gapopen=2, gapextend=2) are more
    lenient than blastn's built-in defaults (2/-3, 5/2) and suit sequences
    with moderate divergence (~95% identity).  See module docstring for the
    full NCBI table of valid reward/penalty -> gap cost combinations.

    -max_hsps 1 and -max_target_seqs 1 ensure a single unambiguous HSP.
    """
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as qf:
        qf.write(query_fasta)
        query_file = qf.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as sf:
        sf.write(subject_fasta)
        subject_file = sf.name

    try:
        cmd = [
            "blastn",
            "-query",           query_file,
            "-subject",         subject_file,
            "-task",            task,
            "-perc_identity",   str(perc_identity),
            "-reward",          str(reward),
            "-penalty",         str(penalty),
            "-gapopen",         str(gapopen),
            "-gapextend",       str(gapextend),
            "-max_hsps",        "1",
            "-max_target_seqs", "1",
            "-outfmt",          "6 qseqid sseqid pident length mismatch gapopen "
                                "qstart qend sstart send evalue bitscore",
            "-dust",            "no",
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

def best_hit(
    blast_output: str,
    query_len: int,
    perc_length: float = 0.0,
) -> tuple[int, int] | None:
    """
    Parse blastn outfmt-6 output and return (sstart, send) for the first
    (highest-scoring) HSP.  sstart > send indicates a minus-strand hit.

    If perc_length > 0, the HSP is only accepted when query coverage
    (qend - qstart + 1) / query_len * 100 >= perc_length.

    Returns None if there are no hits or the hit fails the length filter.
    """
    for line in blast_output.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 12:
            continue
        try:
            qstart = int(cols[6])
            qend   = int(cols[7])
            sstart = int(cols[8])
            send   = int(cols[9])
        except ValueError:
            continue
        if perc_length > 0.0:
            coverage = (abs(qend - qstart) + 1) / query_len * 100
            if coverage < perc_length:
                return None
        return sstart, send
    return None


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
      rel_start < rel_end  ->  forward strand
      rel_start > rel_end  ->  reverse complement (abs_start > abs_end)
    """
    offset    = region_start - 1
    abs_start = rel_start + offset
    abs_end   = rel_end   + offset
    return abs_start, abs_end


# ---------------------------------------------------------------------------
# Per-record worker (must be top-level for multiprocessing pickling)
# ---------------------------------------------------------------------------

def _process_record(
    record_data: tuple[str, str],
    db: str,
    task: str,
    perc_identity: float,
    perc_length: float,
    reward: int,
    penalty: int,
    gapopen: int,
    gapextend: int,
) -> dict:
    """
    Process a single (record_id, seq_str) tuple.
    Returns a dict with keys:
      record_id  : original header ID
      new_id     : updated header ID, or None on failure
      seq_str    : sequence string (unchanged)
      error      : error/warning message string, or None on success
    """
    record_id, seq_str = record_data
    result = {
        "record_id": record_id,
        "new_id":    None,
        "seq_str":   seq_str,
        "error":     None,
    }

    try:
        accession, reg_start, reg_end = parse_header(record_id)
    except ValueError as exc:
        result["error"] = str(exc)
        return result

    try:
        region_fasta = extract_region(db, accession, reg_start, reg_end)
    except RuntimeError as exc:
        result["error"] = str(exc)
        return result

    query_fasta = f">{record_id}\n{seq_str}\n"
    try:
        blast_out = run_blastn(
            query_fasta, region_fasta,
            task=task,
            perc_identity=perc_identity,
            reward=reward,
            penalty=penalty,
            gapopen=gapopen,
            gapextend=gapextend,
        )
    except RuntimeError as exc:
        result["error"] = str(exc)
        return result

    hit = best_hit(blast_out, len(seq_str), perc_length)
    if hit is None:
        result["error"] = (
            "No BLAST hit — try lowering --perc_identity / --perc_length "
            "or using --task blastn-short."
        )
        return result

    rel_s, rel_e = hit
    abs_s, abs_e = relative_to_absolute(reg_start, rel_s, rel_e)
    result["new_id"] = f"{accession}_{abs_s}-{abs_e}"
    return result


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
            "Remap multifasta coordinates by aligning each sequence to its\n"
            "reference genome region via blastn.\n"
            "Input headers must follow the format:  accessionID_start-end\n\n"
            "Scoring defaults (reward=1, penalty=-2, gapopen=2, gapextend=2)\n"
            "are lenient and suited to sequences with moderate divergence\n"
            "(~95%% identity). See module docstring for valid combinations."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # -- I/O -----------------------------------------------------------------
    io = parser.add_argument_group("I/O")
    io.add_argument(
        "-i", "--input",    required=True,
        help="Input multifasta (headers: accessionID_start-end)",
    )
    io.add_argument(
        "-o", "--output",   required=True,
        help="Output multifasta with updated coordinates",
    )
    io.add_argument(
        "-db", "--database", required=True,
        help="Path to local BLAST nucleotide database (blastdbcmd + blastn)",
    )

    # -- BLAST behaviour -----------------------------------------------------
    blast = parser.add_argument_group("BLAST behaviour")
    blast.add_argument(
        "--task", default="blastn",
        choices=["blastn", "blastn-short", "dc-megablast", "megablast"],
        help="blastn task (default: blastn; use blastn-short for sequences <50 bp)",
    )
    blast.add_argument(
        "--perc_identity", type=float, default=80.0,
        help="Minimum %%identity for BLAST hits (default: 80.0)",
    )
    blast.add_argument(
        "--perc_length", type=float, default=0.0,
        help=(
            "Minimum query coverage as a %% of query length (default: 0.0, "
            "disabled). E.g. 90.0 requires the alignment to span at least "
            "90%% of the query sequence."
        ),
    )

    # -- Scoring -------------------------------------------------------------
    scoring = parser.add_argument_group(
        "Scoring",
        "Reward/penalty ratio controls mismatch tolerance: 1/-1 ~75%% id, "
        "1/-2 ~95%% id, 2/-3 ~99%% id (blastn default). "
        "Gap costs must form a valid combination with the chosen ratio "
        "(see module docstring).",
    )
    scoring.add_argument(
        "--reward",    type=int, default=1,
        help="Match reward (default: 1)",
    )
    scoring.add_argument(
        "--penalty",   type=int, default=-2,
        help="Mismatch penalty, must be negative (default: -2)",
    )
    scoring.add_argument(
        "--gapopen",   type=int, default=2,
        help="Gap opening cost (default: 2)",
    )
    scoring.add_argument(
        "--gapextend", type=int, default=2,
        help="Gap extension cost (default: 2)",
    )

    # -- Runtime -------------------------------------------------------------
    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of parallel worker processes (default: 1)",
    )
    runtime.add_argument(
        "--keep_on_fail", action="store_true",
        help=(
            "Retain records unchanged in output when alignment fails. "
            "Default: omit failed records and report to stderr."
        ),
    )

    args = parser.parse_args()

    # -- Load all records upfront --------------------------------------------
    with open(args.input) as fh:
        all_records = list(SeqIO.parse(fh, "fasta"))

    total = len(all_records)
    print(
        f"[INFO]  {total} records loaded. Running with {args.threads} thread(s).\n"
        f"[INFO]  Scoring: reward={args.reward}, penalty={args.penalty}, "
        f"gapopen={args.gapopen}, gapextend={args.gapextend}",
        file=sys.stderr,
    )

    record_data = [(rec.id, str(rec.seq)) for rec in all_records]

    # Bind fixed arguments to the worker using functools.partial.
    # _process_record must be a top-level function for pickle compatibility.
    worker = partial(
        _process_record,
        db=args.database,
        task=args.task,
        perc_identity=args.perc_identity,
        perc_length=args.perc_length,
        reward=args.reward,
        penalty=args.penalty,
        gapopen=args.gapopen,
        gapextend=args.gapextend,
    )

    # -- Run in parallel (Pool(1) is valid and incurs negligible overhead) ---
    with Pool(processes=args.threads) as pool:
        results = list(pool.imap(worker, record_data))

    # -- Collect output and report -------------------------------------------
    records_out: list[SeqRecord] = []
    n_updated = n_failed = 0

    for res in results:
        if res["error"] is None:
            print(f"[INFO]  {res['record_id']}  ->  {res['new_id']}", file=sys.stderr)
            records_out.append(
                SeqRecord(Seq(res["seq_str"]), id=res["new_id"], description="", name="")
            )
            n_updated += 1
        else:
            print(f"[WARNING]  {res['record_id']}: {res['error']}", file=sys.stderr)
            n_failed += 1
            if args.keep_on_fail:
                records_out.append(
                    SeqRecord(
                        Seq(res["seq_str"]),
                        id=res["record_id"],
                        description="",
                        name="",
                    )
                )

    # -- Write output --------------------------------------------------------
    write_fasta(records_out, args.output)
    print(
        f"\n[DONE]  {n_updated} records updated, {n_failed} failed. "
        f"Output -> {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
