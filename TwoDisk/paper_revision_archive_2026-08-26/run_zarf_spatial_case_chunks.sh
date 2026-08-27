#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 5 ]]; then
    echo "Usage: $0 MASS_FLUX [CANDIDATE_COUNT] [WORKER_INDEX] [WORKER_COUNT] [THREADS]" >&2
    exit 2
fi

mass_flux="$1"
candidate_count="${2:-1}"
worker_index="${3:-0}"
worker_count="${4:-1}"
threads="${5:-10}"
((worker_index >= 0 && worker_index < worker_count)) || {
    echo "WORKER_INDEX must satisfy 0 <= index < WORKER_COUNT" >&2
    exit 2
}
root="/home/zhj/Rotating-Flow-ToolKit/TwoDisk"
zarf_root="/home/zhj/main/code/compress/compressible_bek/TwoDisk/data/stability/Zarf_neutral"
julia="/home/zhj/.juliaup/bin/julia"
temporal_directory=$(printf "%s/Ts=%.1f" "$zarf_root" "$mass_flux")
case_label=$(printf "%+.1f" "$mass_flux" | sed -e 's/+//g' -e 's/-/m/g' -e 's/\./p/g')
if [[ "$mass_flux" != -* ]]; then
    case_label="p${case_label}"
fi
case_directory="Ts_${case_label}_N99"
output_root="$root/zarf_spatial_growth_results/remaining_newton_chunks_c${candidate_count}/Ts_${case_label}"
mkdir -p "$output_root"

mapfile -t radii < <(
    find "$temporal_directory" -maxdepth 1 -type f -name 'R=*.dat' -printf '%f\n' |
    sed -e 's/^R=//' -e 's/\.dat$//' |
    sort -g
)

chunk_index=0
for ((start = 0; start < ${#radii[@]}; start += 10)); do
    if ((chunk_index % worker_count != worker_index)); then
        ((chunk_index += 1))
        continue
    fi
    end=$((start + 9))
    ((end >= ${#radii[@]})) && end=$((${#radii[@]} - 1))
    radius_min="${radii[start]}"
    radius_max="${radii[end]}"
    chunk_label=$(printf 'R%s_%s' "$radius_min" "$radius_max" | tr '.' 'p')
    chunk_directory="$output_root/$chunk_label"
    envelope="$chunk_directory/$case_directory/spatial_growth_envelope.dat"
    if [[ -s "$envelope" ]]; then
        echo "SKIP a_s=$mass_flux $chunk_label"
        ((chunk_index += 1))
        continue
    fi
    echo "START a_s=$mass_flux $chunk_label $(date '+%F %T')"
    "$julia" -t "$threads" --project="$root" "$root/zarf_spatial_growth_scan.jl" \
        --ts "$mass_flux" \
        --R-min "$radius_min" \
        --R-max "$radius_max" \
        --frequency-samples 13 \
        --candidate-count "$candidate_count" \
        --output-root "$chunk_directory" \
        2>&1 | tee "$output_root/$chunk_label.log"
    echo "DONE a_s=$mass_flux $chunk_label $(date '+%F %T')"
    ((chunk_index += 1))
done
