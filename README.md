# ANI Calculator

ANI Calculator is a Docker-first wrapper around FastANI for computing pairwise whole-genome Average Nucleotide Identity (ANI) across multiple nucleotide FASTA files in a folder.

The pipeline scans top-level FASTA files, runs all unique genome pairs, and writes:
- a long-format pairwise ANI table
- a square ANI similarity matrix

## Original Tool and Container

This repository uses the following upstream resources:
- FastANI (original ANI engine): https://github.com/ParBLiSS/FastANI
- Default Docker image used by this pipeline: `staphb/fastani:1.34`
- Docker image page: https://hub.docker.com/r/staphb/fastani

## Requirements

- Docker installed and daemon running
- Bash shell
- Input genomes as nucleotide FASTA files

## Quick Usage

Run with defaults (auto-loads `ani_reporting.conf` if present next to the script):

```bash
./ANI_script.sh -i dna_sequences -o results/my_run
```

Run with an explicit config file:

```bash
./ANI_script.sh -i dna_sequences -o results/my_run -c /path/to/ani_reporting.conf
```

Optional Docker image override:

```bash
./ANI_script.sh -i dna_sequences -o results/my_run -d staphb/fastani:1.34
```

## Input and Output

Input expectations:
- Top-level files only (no recursive scan)
- Accepted extensions: `.fa`, `.fna`, `.fasta`
- Nucleotide alphabet expected: A, T, C, G, N
- At least 2 FASTA files are required
- Filename stems must be unique (used as species IDs)

Outputs for `-o results/my_run`:
- `results/my_run_long.tsv`: pairwise rows with `species_a`, `species_b`, `ANI`, `fragments_mapped`, `total_fragments`
- `results/my_run_matrix.tsv`: square species-by-species ANI matrix

## Configuration (`ani_reporting.conf`)

The config file controls ANI reporting defaults without editing `ANI_script.sh`.

Supported keys:
- `DOCKER_IMAGE`: container image used to run FastANI
- `THREADS`: FastANI thread count (`-t`)
- `MIN_ANI`: minimum ANI required for numeric reporting
- `MIN_ALIGNMENT_FRACTION`: minimum mapping fraction required, computed as `fragments_mapped / total_fragments`
- `MIN_MAPPED_FRAGMENTS`: minimum mapped fragment count required
- `BELOW_THRESHOLD_ACTION`: behavior for pairs failing thresholds (`keep`, `na`, `drop`)
- `NA_VALUE`: placeholder value written for filtered or missing ANI values

Behavior notes:
- `keep`: keep FastANI ANI value even if thresholds fail
- `na`: write `NA_VALUE` for ANI when thresholds fail
- `drop`: omit failed pairs from long table (matrix cells remain `NA_VALUE` if no accepted pair exists)

## Why These Default Values

Current defaults in `ani_reporting.conf` are conservative and intended for reproducible species-level screening:
- `MIN_ANI=80.0`: practical floor for reporting species-similarity ANI in many workflows
- `MIN_ALIGNMENT_FRACTION=0.20`: requires at least modest query coverage before trusting ANI
- `MIN_MAPPED_FRAGMENTS=50`: avoids reporting from very sparse fragment support
- `BELOW_THRESHOLD_ACTION=na`: keeps pair structure in outputs while clearly flagging low-confidence pairs
- `DOCKER_IMAGE=staphb/fastani:1.34`: pinned default image for reproducibility across runs

These are defaults, not hard biological rules. Adjust for your dataset quality and analysis goals.

## Notes and Limitations

- FastANI can report asymmetric values depending on query/reference direction.
- FastANI may still report some ANI values below 80 in practice; this wrapper applies post-run threshold filtering based on config.
- Low assembly quality can increase missing or filtered results.
