#!/usr/bin/env bash
# Build a publication-ready bacterial genome similarity tree from whole-genome
# FASTA/FNA files using Mashtree, then reroot with newick-utils.
#
# This script is aimed at practical visualization in proteomics/MS workflows,
# not deep evolutionary reconstruction.

set -euo pipefail
IFS=$'\n\t'

usage() {
    cat <<'EOF'
Usage:
  ./mashtree_phylo_pipeline.sh --install
  ./mashtree_phylo_pipeline.sh -i INPUT_DIR [-t THREADS] [-o OUTPUT_PREFIX]

Options:
  --install           Install required Ubuntu dependencies and tools
  -i INPUT_DIR        Directory containing .fna/.fasta genome files (top-level)
  -t THREADS          Thread count for mashtree --numcpus (default: nproc)
  -o OUTPUT_PREFIX    Output prefix (default: genomes_tree)
  -h, --help          Show this help

Outputs:
  <prefix>.nwk
  <prefix>_rooted.nwk

Example:
  ./mashtree_phylo_pipeline.sh -i "dna_sequences" -t 16 -o genomes_tree
EOF
}

log() { echo "INFO: $*" >&2; }
err() { echo "ERROR: $*" >&2; }
have_cmd() { command -v "$1" >/dev/null 2>&1; }

install_newick_utils_from_source() {
    log "Installing newick-utils from source"
    sudo apt-get install -y build-essential autoconf automake libtool flex bison git libxml2-dev

    local tmp_dir
    tmp_dir="$(mktemp -d)"
    trap 'rm -rf "${tmp_dir}"' RETURN

    git clone --depth 1 https://github.com/tjunier/newick_utils "${tmp_dir}/newick_utils"
    cd "${tmp_dir}/newick_utils"

    # Avoid inherited Conda/Miniforge include or linker flags and force
    # -fcommon for older newick-utils sources on modern GCC toolchains.
    env -u CPATH -u C_INCLUDE_PATH -u CPLUS_INCLUDE_PATH \
        -u LIBRARY_PATH -u LD_LIBRARY_PATH -u CPPFLAGS -u LDFLAGS \
        autoreconf -fi
    env -u CPATH -u C_INCLUDE_PATH -u CPLUS_INCLUDE_PATH \
        -u LIBRARY_PATH -u LD_LIBRARY_PATH \
        CPPFLAGS= LDFLAGS= CFLAGS="-O2 -fcommon" \
        ./configure --prefix=/usr/local
    env -u CPATH -u C_INCLUDE_PATH -u CPLUS_INCLUDE_PATH \
        -u LIBRARY_PATH -u LD_LIBRARY_PATH \
        CPPFLAGS= LDFLAGS= CFLAGS="-O2 -fcommon" \
        make -j"$(nproc)"
    sudo make install
}

install_mashtree_from_source() {
    log "Installing Mashtree via cpanm"
    sudo apt-get install -y perl cpanminus mash quicktree git
    sudo cpanm --notest Bio::Perl
    sudo cpanm --notest https://github.com/lskatz/mashtree/archive/refs/heads/master.tar.gz
}

install_dependencies() {
    log "Installing Ubuntu dependencies"
    sudo apt-get update
    sudo apt-get install -y bash findutils coreutils mawk perl cpanminus mash quicktree git

    if ! have_cmd mashtree; then
        install_mashtree_from_source
    fi

    if ! have_cmd nw_reroot; then
        if ! sudo apt-get install -y newick-utils; then
            install_newick_utils_from_source
        fi
    fi

    have_cmd mashtree || { err "mashtree not found after install"; exit 1; }
    have_cmd mash || { err "mash not found after install"; exit 1; }
    have_cmd nw_reroot || { err "nw_reroot not found after install"; exit 1; }
    log "Dependency installation completed"
}

install_mode=0
input_dir=""
threads="$(nproc)"
output_prefix="genomes_tree"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --install)
            install_mode=1
            shift
            ;;
        -i)
            input_dir="${2:-}"
            shift 2
            ;;
        -t)
            threads="${2:-}"
            shift 2
            ;;
        -o)
            output_prefix="${2:-}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            err "Unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ "${install_mode}" -eq 1 ]]; then
    install_dependencies
    exit 0
fi

[[ -n "${input_dir}" ]] || { err "-i INPUT_DIR is required"; usage; exit 1; }
[[ -d "${input_dir}" ]] || { err "Input directory not found: ${input_dir}"; exit 1; }
[[ "${threads}" =~ ^[1-9][0-9]*$ ]] || { err "Threads must be a positive integer"; exit 1; }

have_cmd mashtree || { err "mashtree not found. Run with --install first."; exit 1; }
have_cmd nw_reroot || { err "nw_reroot not found. Run with --install first."; exit 1; }

# Discover top-level FASTA/FNA inputs safely with null delimiters.
declare -a genomes=()
while IFS= read -r -d '' file; do
    genomes+=("${file}")
done < <(find "${input_dir}" -maxdepth 1 -type f \( -iname '*.fna' -o -iname '*.fasta' \) -print0 | sort -z)

count="${#genomes[@]}"
if [[ "${count}" -lt 2 ]]; then
    err "Need at least 2 .fna/.fasta files in ${input_dir}. Found: ${count}"
    exit 1
fi

out_dir="$(dirname "${output_prefix}")"
mkdir -p "${out_dir}"

unrooted_tree="${output_prefix}.nwk"
rooted_tree="${output_prefix}_rooted.nwk"

log "Found ${count} genome files"
log "Running Mashtree with --numcpus ${threads}"

# Generate unrooted Newick tree from whole-genome assemblies.
mashtree --numcpus "${threads}" "${genomes[@]}" > "${unrooted_tree}"

[[ -s "${unrooted_tree}" ]] || { err "Mashtree produced empty output: ${unrooted_tree}"; exit 1; }
log "Unrooted tree: ${unrooted_tree}"

# Root with newick-utils. With no outgroup label, nw_reroot reroots on the
# longest branch; this is used here as a simple midpoint-style publication view.
nw_reroot "${unrooted_tree}" > "${rooted_tree}"

[[ -s "${rooted_tree}" ]] || { err "Rooted tree output is empty: ${rooted_tree}"; exit 1; }
log "Rooted tree: ${rooted_tree}"
log "Done. Upload rooted tree to iTOL or FigTree."
