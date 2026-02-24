## **Design Choices**

## **Why BLAST over pairwise alignment?**

Your sequences have mismatches and could be on either strand — BLAST handles gapped, mismatch-tolerant local alignment natively. Since you already have a local BLAST db, using `blastdbcmd` + `blastn -subject` keeps everything self-contained without needing to load full genome sequences into memory.

## **Workflow per record**

1.  **Parse header** — `rsplit('_', 1)` splits on the *last* underscore only, so accession IDs like `NZ_CP012345.1_100-200` are handled correctly

2.  **Extract reference region** — `blastdbcmd -range start-end` pulls only the relevant window from the genome (1-based, inclusive)​

3.  **Run `blastn -subject`** — aligns query against the extracted region using temp files; `--dust no` avoids masking IS-element-relevant repetitive sequence

4.  **Parse best HSP** — selects the hit with the highest bitscore from `-outfmt 6`

5.  **Remap coordinates** — relative BLAST subject coords → absolute genome coords via: `abs = rel + (region_start - 1)`

## **Coordinate conversion**

For a region extracted starting at genome position RR, and a BLAST `sstart`/`send` of pp:

abs=p+(R−1)abs=p+(R−1)

Strand orientation is preserved: if `sstart > send` in the output header, the sequence hit on the minus strand.

------------------------------------------------------------------------

## **Usage**

```         
bash
```

`# Standard use
python remap_fasta_coords.py \
    -i input.fasta \
    -o output_remapped.fasta \
    -db /path/to/blast/db

# Short sequences (<50 bp, e.g. small IS target sites)
python remap_fasta_coords.py \
    -i input.fasta \
    -o output_remapped.fasta \
    -db /path/to/blast/db \
    --task blastn-short \
    --perc_identity 70

# Keep failed records in output unchanged (rather than dropping them)
python remap_fasta_coords.py -i input.fasta -o out.fasta -db mydb --keep_on_fail`

------------------------------------------------------------------------

## **Key Arguments**

| **Flag** | **Default** | **Purpose** |
|:---|:---|:---|
| `-i` / `--input` | required | Input multifasta |
| `-o` / `--output` | required | Output multifasta |
| `-db` / `--database` | required | BLAST db path (used by both `blastdbcmd` and `blastn`) |
| `--task` | `blastn` | Switch to `blastn-short` for sequences under \~50 bp |
| `--perc_identity` | `80.0` | Lower if sequences are highly diverged |
| `--keep_on_fail` | off | Retain original record in output on failure instead of omitting |

------------------------------------------------------------------------

## **Dependencies**

-   **BLAST+** (`blastn`, `blastdbcmd`) — must be on `$PATH`

-   **Biopython** — `pip install biopython` or `conda install biopython`

One thing to note: if your IS element sequences are on the minus strand relative to the reference, the output header will have `abs_start > abs_end` (e.g. `NZ_CP012345.1_2500-2200`). This is intentional and consistent with how tools like NCBI and many annotation pipelines encode strand. If you'd prefer separate `+/-` strand notation instead, let me know and I can add a `--strand_notation` flag.

Documentation and code logic was prepared with the assistance of Claude Sonnet 4.6 Thinking via Perplexity Pro