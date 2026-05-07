#!/usr/bin/env bash
# This script utilizes the FastANI tool hosted at: https://github.com/ParBLiSS/FastANI
#
# Usage: ANI_script.sh -i <input_dir> -o <output_prefix> [-d <docker_image>] [-c <config_file>]
#
# Input FASTA expectations:
#   - Nucleotide sequences only (A, T, C, G, N accepted; full IUPAC out of scope)
#   - N bases are accepted silently; handling delegated to FastANI
#   - One assembled genome per file (species-level comparison)
#   - Top-level files only (.fa, .fna, .fasta); no recursion, no compressed inputs
#   - Minimum 2 FASTA files required
#   - Filename stem used as species identifier; stems must be unique
#
# Outputs:
#   <output_prefix>_long.tsv    — pairwise ANI long table
#   <output_prefix>_matrix.tsv  — square species × species ANI matrix

set -uo pipefail

# ── Constants ─────────────────────────────────────────────────────────────────
readonly DEFAULT_IMAGE="staphb/fastani:1.34"
readonly DEFAULT_CONFIG_FILE="ani_reporting.conf"

# Runtime defaults (can be overridden by config)
DOCKER_IMAGE="${DEFAULT_IMAGE}"
THREADS="1"
MIN_ANI="80.0"
MIN_ALIGNMENT_FRACTION="0.20"
MIN_MAPPED_FRAGMENTS="50"
BELOW_THRESHOLD_ACTION="na"
NA_VALUE="NA"

float_lt() {
    awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 < b + 0) }'
}

validate_decimal() {
    [[ "$1" =~ ^[0-9]+([.][0-9]+)?$ ]]
}

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat >&2 <<EOF
Usage: $(basename "$0") -i <input_dir> -o <output_prefix> [-d <docker_image>] [-c <config_file>]

  -i <input_dir>      Directory containing nucleotide FASTA files (.fa, .fna, .fasta)
  -o <output_prefix>  Prefix for output files (e.g. results/my_run)
  -d <docker_image>   FastANI Docker image (default: ${DEFAULT_IMAGE})
    -c <config_file>    Optional config file for ANI reporting thresholds
  -h                  Show this help message

Config keys supported:
    MIN_ANI, MIN_ALIGNMENT_FRACTION, MIN_MAPPED_FRAGMENTS,
    BELOW_THRESHOLD_ACTION (keep|na|drop), NA_VALUE, DOCKER_IMAGE, THREADS
EOF
    exit 1
}

# ── Argument parsing ──────────────────────────────────────────────────────────
input_dir=""
output_prefix=""
docker_image_cli=""
config_file=""

while getopts ":i:o:d:c:h" opt; do
    case "${opt}" in
        i) input_dir="${OPTARG}" ;;
        o) output_prefix="${OPTARG}" ;;
        d) docker_image_cli="${OPTARG}" ;;
        c) config_file="${OPTARG}" ;;
        h) usage ;;
        :) echo "ERROR: Option -${OPTARG} requires an argument." >&2; usage ;;
        \?) echo "ERROR: Unknown option -${OPTARG}." >&2; usage ;;
    esac
done

[[ -z "${input_dir}" ]]     && { echo "ERROR: -i <input_dir> is required." >&2;     usage; }
[[ -z "${output_prefix}" ]] && { echo "ERROR: -o <output_prefix> is required." >&2; usage; }

# Strip trailing slashes so -o results/ and -o results behave identically
output_prefix="${output_prefix%/}"

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
default_config_path="${script_dir}/${DEFAULT_CONFIG_FILE}"

if [[ -z "${config_file}" && -f "${default_config_path}" ]]; then
    config_file="${default_config_path}"
fi

if [[ -n "${config_file}" ]]; then
    if [[ ! -f "${config_file}" ]]; then
        echo "ERROR: Config file '${config_file}' not found." >&2
        exit 1
    fi
    # shellcheck source=/dev/null
    source "${config_file}"
    echo "INFO: Loaded config file '${config_file}'."
fi

if [[ -n "${docker_image_cli}" ]]; then
    docker_image="${docker_image_cli}"
else
    docker_image="${DOCKER_IMAGE}"
fi

if ! validate_decimal "${MIN_ANI}"; then
    echo "ERROR: MIN_ANI must be a decimal number. Got '${MIN_ANI}'." >&2
    exit 1
fi
if ! awk -v v="${MIN_ANI}" 'BEGIN { exit !(v >= 0 && v <= 100) }'; then
    echo "ERROR: MIN_ANI must be between 0 and 100. Got '${MIN_ANI}'." >&2
    exit 1
fi

if ! validate_decimal "${MIN_ALIGNMENT_FRACTION}"; then
    echo "ERROR: MIN_ALIGNMENT_FRACTION must be a decimal number. Got '${MIN_ALIGNMENT_FRACTION}'." >&2
    exit 1
fi
if ! awk -v v="${MIN_ALIGNMENT_FRACTION}" 'BEGIN { exit !(v >= 0 && v <= 1) }'; then
    echo "ERROR: MIN_ALIGNMENT_FRACTION must be between 0 and 1. Got '${MIN_ALIGNMENT_FRACTION}'." >&2
    exit 1
fi

if [[ ! "${MIN_MAPPED_FRAGMENTS}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: MIN_MAPPED_FRAGMENTS must be a non-negative integer. Got '${MIN_MAPPED_FRAGMENTS}'." >&2
    exit 1
fi

if [[ ! "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: THREADS must be a positive integer. Got '${THREADS}'." >&2
    exit 1
fi

case "${BELOW_THRESHOLD_ACTION}" in
    keep|na|drop) ;;
    *)
        echo "ERROR: BELOW_THRESHOLD_ACTION must be one of keep, na, drop. Got '${BELOW_THRESHOLD_ACTION}'." >&2
        exit 1
        ;;
esac

[[ ! -d "${input_dir}" ]]   && { echo "ERROR: Input directory '${input_dir}' not found." >&2; exit 1; }

# ── Dependency checks ─────────────────────────────────────────────────────────
if ! command -v docker &>/dev/null; then
    echo "ERROR: Docker is not installed or not in PATH." >&2
    exit 1
fi

if ! docker info &>/dev/null 2>&1; then
    echo "ERROR: Docker daemon is not running or current user lacks permission." >&2
    exit 1
fi

# ── Input discovery ───────────────────────────────────────────────────────────
mapfile -t fasta_files < <(
    find "${input_dir}" -maxdepth 1 -type f \
        \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | sort
)

n_files="${#fasta_files[@]}"

if [[ "${n_files}" -lt 2 ]]; then
    echo "ERROR: At least 2 FASTA files are required in '${input_dir}'. Found: ${n_files}." >&2
    exit 1
fi

# ── Validate files and derive species names ───────────────────────────────────
declare -A stem_seen
species_names=()

for f in "${fasta_files[@]}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: File '${f}' is empty or zero-byte." >&2
        exit 1
    fi

    first_char=$(head -c 1 "${f}")
    if [[ "${first_char}" != ">" ]]; then
        echo "ERROR: File '${f}' does not appear to be a valid FASTA (first character is not '>')." >&2
        exit 1
    fi

    base=$(basename "${f}")
    stem="${base%.*}"

    if [[ -n "${stem_seen[${stem}]+_}" ]]; then
        echo "ERROR: Duplicate filename stem '${stem}' detected. Species names must be unique." >&2
        exit 1
    fi

    stem_seen["${stem}"]=1
    species_names+=("${stem}")
done

echo "INFO: Found ${n_files} FASTA files."
echo "INFO: Species: ${species_names[*]}"
echo "INFO: Using Docker image: ${docker_image}"
echo "INFO: Reporting thresholds: MIN_ANI=${MIN_ANI}, MIN_ALIGNMENT_FRACTION=${MIN_ALIGNMENT_FRACTION}, MIN_MAPPED_FRAGMENTS=${MIN_MAPPED_FRAGMENTS}, BELOW_THRESHOLD_ACTION=${BELOW_THRESHOLD_ACTION}"

# ── Temp workspace ────────────────────────────────────────────────────────────
tmp_dir=$(mktemp -d)
trap 'rm -rf "${tmp_dir}"' EXIT

input_abs=$(realpath "${input_dir}")

# Ensure output directory exists
output_dir=$(dirname "${output_prefix}")
mkdir -p "${output_dir}"

# ── Pairwise ANI computation ──────────────────────────────────────────────────
long_table="${output_prefix}_long.tsv"
printf "species_a\tspecies_b\tANI\tfragments_mapped\ttotal_fragments\n" > "${long_table}"

n_pairs=$(( n_files * (n_files - 1) / 2 ))
pair_num=0
reported_pairs=0
threshold_na_pairs=0
dropped_pairs=0
no_result_pairs=0

for (( i=0; i<n_files; i++ )); do
    for (( j=i+1; j<n_files; j++ )); do
        query_file="${fasta_files[$i]}"
        ref_file="${fasta_files[$j]}"
        query_stem="${species_names[$i]}"
        ref_stem="${species_names[$j]}"
        pair_num=$(( pair_num + 1 ))

        echo "INFO: [${pair_num}/${n_pairs}] ${query_stem} vs ${ref_stem}"

        query_name=$(basename "${query_file}")
        ref_name=$(basename "${ref_file}")
        out_name="${query_stem}__${ref_stem}.txt"

        if docker run --rm \
                -v "${input_abs}:/data:ro" \
                -v "${tmp_dir}:/out" \
                "${docker_image}" \
                fastANI \
                    -q "/data/${query_name}" \
                    -r "/data/${ref_name}" \
                    -t "${THREADS}" \
                    -o "/out/${out_name}" \
                2>/dev/null; then

            result_file="${tmp_dir}/${out_name}"
            if [[ -f "${result_file}" && -s "${result_file}" ]]; then
                read -r _ _ ani frags total < "${result_file}"

                alignment_fraction=$(awk -v m="${frags}" -v t="${total}" 'BEGIN { if (t == 0) print 0; else print m / t }')
                fails_threshold=0

                if float_lt "${ani}" "${MIN_ANI}"; then
                    fails_threshold=1
                fi
                if float_lt "${alignment_fraction}" "${MIN_ALIGNMENT_FRACTION}"; then
                    fails_threshold=1
                fi
                if (( frags < MIN_MAPPED_FRAGMENTS )); then
                    fails_threshold=1
                fi

                if (( fails_threshold == 0 )); then
                    printf "%s\t%s\t%s\t%s\t%s\n" \
                        "${query_stem}" "${ref_stem}" "${ani}" "${frags}" "${total}" >> "${long_table}"
                    reported_pairs=$(( reported_pairs + 1 ))
                else
                    case "${BELOW_THRESHOLD_ACTION}" in
                        keep)
                            printf "%s\t%s\t%s\t%s\t%s\n" \
                                "${query_stem}" "${ref_stem}" "${ani}" "${frags}" "${total}" >> "${long_table}"
                            reported_pairs=$(( reported_pairs + 1 ))
                            ;;
                        na)
                            printf "%s\t%s\t%s\t%s\t%s\n" \
                                "${query_stem}" "${ref_stem}" "${NA_VALUE}" "${frags}" "${total}" >> "${long_table}"
                            threshold_na_pairs=$(( threshold_na_pairs + 1 ))
                            ;;
                        drop)
                            dropped_pairs=$(( dropped_pairs + 1 ))
                            ;;
                    esac
                fi
            else
                # FastANI produced no output — similarity below threshold
                if [[ "${BELOW_THRESHOLD_ACTION}" == "drop" ]]; then
                    dropped_pairs=$(( dropped_pairs + 1 ))
                else
                    printf "%s\t%s\t%s\t%s\t%s\n" "${query_stem}" "${ref_stem}" "${NA_VALUE}" "${NA_VALUE}" "${NA_VALUE}" >> "${long_table}"
                    threshold_na_pairs=$(( threshold_na_pairs + 1 ))
                fi
                no_result_pairs=$(( no_result_pairs + 1 ))
                echo "WARNING: No ANI result for ${query_stem} vs ${ref_stem} (below similarity threshold)."
            fi
        else
            if [[ "${BELOW_THRESHOLD_ACTION}" == "drop" ]]; then
                dropped_pairs=$(( dropped_pairs + 1 ))
            else
                printf "%s\t%s\t%s\t%s\t%s\n" "${query_stem}" "${ref_stem}" "${NA_VALUE}" "${NA_VALUE}" "${NA_VALUE}" >> "${long_table}"
                threshold_na_pairs=$(( threshold_na_pairs + 1 ))
            fi
            echo "WARNING: Docker command failed for ${query_stem} vs ${ref_stem}. Recorded as NA."
        fi
    done
done

echo "INFO: Long table written to '${long_table}'."

# ── Square matrix generation ──────────────────────────────────────────────────
matrix_file="${output_prefix}_matrix.tsv"

# Build pipe-delimited species list for awk (safe against spaces in stems)
species_list=$(IFS="|"; echo "${species_names[*]}")

awk -v OFS='\t' -v species_list="${species_list}" -v na_value="${NA_VALUE}" '
BEGIN {
    n = split(species_list, sp, "|")
    for (i = 1; i <= n; i++) species[i] = sp[i]
}
NR > 1 && $3 != na_value {
    key[$1 SUBSEP $2] = $3
    key[$2 SUBSEP $1] = $3
}
END {
    # Header row
    printf ""
    for (i = 1; i <= n; i++) printf "\t%s", species[i]
    printf "\n"
    # Data rows
    for (i = 1; i <= n; i++) {
        printf "%s", species[i]
        for (j = 1; j <= n; j++) {
            if (i == j) {
                printf "\t100"
            } else if ((species[i] SUBSEP species[j]) in key) {
                printf "\t%s", key[species[i] SUBSEP species[j]]
            } else {
                printf "\t%s", na_value
            }
        }
        printf "\n"
    }
}
' "${long_table}" > "${matrix_file}"

echo "INFO: Square matrix written to '${matrix_file}'."

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " ANI pipeline complete"
echo " Species:        ${n_files}"
echo " Pairs computed: ${n_pairs}"
echo " Numeric ANI:    ${reported_pairs}"
echo " Threshold -> NA:${threshold_na_pairs}"
echo " Dropped pairs:  ${dropped_pairs}"
echo " No-result pairs:${no_result_pairs}"
echo " Long table:     ${long_table}"
echo " Matrix:         ${matrix_file}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

