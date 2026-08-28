# CAUSAL-flow

Nextflow DSL2 pipeline for colocalization (`coloc`) and Mendelian Randomization analysis between
GWAS and eQTL datasets.

## Rules

1. Ask, don't assume. If something is unclear, ask before writing a single line.
2. Simplest solution first. No abstractions or flexibility beyond what was asked.
3. Don't touch unrelated code. If it's not part of the task, leave it alone.
4. Flag uncertainty explicitly. If not confident, say so before proceeding.
5. Before any significant task, show 2-3 approaches and wait for a choice.
6. Never open with filler phrases. Start with the answer.
7. Match response length to task complexity. No padding, no restatements.

## Commits

Never add a "Co-Authored-By: Claude" / "Claude-Session" trailer, or any other Claude/Anthropic
attribution, to commits in this repo — Alex is the sole author of record on every commit, past and
future. If a commit template or tool default would normally add one, omit it here.

## Memory and decisions

Maintain `MEMORY.md` in `.claude/`. After any significant decision, add a dated entry: what was
decided / why / what was rejected and why. Read `MEMORY.md` at the start of every session.

## Who I am

Alex — Computational Biologist. Strong in: single-cell, QTL, coloc, MR, bioinformatics pipelines.
Still learning: deep single-cell analysis, AI/ML.

Do not over-explain concepts I know. Assume familiarity with GWAS summary stats, eQTL mapping,
LD clumping, Bayesian colocalization (coloc), Mendelian Randomization, and Nextflow DSL2.

## Default behaviors

- No comments unless the WHY is non-obvious.
- No new files unless necessary — edit existing ones.

## Architecture

- `main.nf` — entry point; wires GWAS/eQTL channel construction (single file or directory + glob
  pattern) through the pipeline steps below.
- `modules/clump_gwas.nf` — LD-based clumping of GWAS summary stats (wraps `R/clump_gwas_source.R`).
- `modules/run_coloc.nf` — `run_coloc` (per GWAS × eQTL pair colocalization) + `combine_coloc`
  (aggregates results, filters on `--coloc_pph4_threshold`). Wraps `R/run_coloc_functions.R`.
- `modules/run_MR.nf` — optional MR step: `run_MR` + `combine_MR`, only runs when
  `--allele_file` is supplied. Wraps `R/run_MR_functions.R`.
- `R/` — the actual analysis code the modules call into (`clump_gwas_source.R`,
  `run_coloc_functions.R`, `run_MR_functions.R`).
- `dockerfile`, `build_docker.sh` — container image (`ghcr.io/johnsonlab-ic/causal-flow`).
- `nextflow.config` — profiles: `standard` (local, Docker), `offline` (local, Docker, no network),
  `imperial` (PBS Pro + Singularity, resource labels `process_single/low/medium/high/long/high_memory`).

## Pipeline flow

```
GWAS file(s) → CLUMP_GWAS (LD clumping)
                                          ⇘
eQTL file(s) ───────────────────────────── combine → RUN_COLOC (per pair) → COMBINE_COLOC (PP.H4 filter)
                                                                                  ⇘ [if --allele_file]
                                                                             RUN_MR → COMBINE_MR
```

### Key params

- `--gwas_file` / `--gwas_dir` + `--gwas_pattern` — single file or directory of GWAS summary stats
- `--eqtl_file` / `--eqtl_dir` + `--eqtl_pattern` — single file or directory of eQTL data
- `--gene_location_file` — required, gene coordinates for coloc window construction
- `--pvalue` (default `5e-8`), `--window_size` (default `1e6`) — clumping thresholds
- `--coloc_pph4_threshold` (default `0.2`) — PP.H4 filter for combined coloc results
- `--allele_file` — presence gates the optional MR step
- `--eqtl_fdr_threshold` (default `0.05`), `--mr_pval_threshold` (default `0.05`) — MR filters
- `--outdir` — output directory

## Known gaps

- `manifest.author` in `nextflow.config` still has a placeholder email — fix when touching that file.
- `PLOT_RESULTS` / `GENERATE_REPORT` modules are referenced (commented out) in `main.nf` but not
  yet implemented.
