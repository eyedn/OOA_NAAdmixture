#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           map_slurm_task.sh
###############################################################################

# shared Slurm task-index mappers for chromosome, replicate, and population.

##### main functions ##########################################################
# map a chromosome-major Slurm task ID to chromosome and replicate
map_slurm_task_to_chr_rep() {
    local task_id="$1"
    local num_reps="$2"
    shift 2
    if [[ "${1:-}" != "--" ]]; then
        echo "ERROR: expected -- before chromosome array" >&2
        return 1
    fi
    shift
    local chroms=( "$@" )
    local num_chroms="${#chroms[@]}"
    local max_task=$((num_reps * num_chroms))
    local chr_index
    local rep
    local chr

    if (( num_reps < 1 )); then
        echo "ERROR: num_reps must be >= 1" >&2
        return 1
    fi
    if (( num_chroms < 1 )); then
        echo "ERROR: at least one chromosome is required" >&2
        return 1
    fi
    if (( task_id < 1 || task_id > max_task )); then
        echo \
            "ERROR: task ID ${task_id} outside valid range 1..${max_task}" \
            >&2
        return 1
    fi

    # each chromosome owns a contiguous, replicate-minor task block.
    chr_index=$(((task_id - 1) / num_reps))
    rep=$((((task_id - 1) % num_reps) + 1))
    chr="${chroms[chr_index]}"

    printf "%s %s\n" "${chr}" "${rep}"
}


# map a chromosome-major Slurm task ID to chromosome and population
map_slurm_task_to_chr_pop() {
    local task_id="$1"
    shift
    if [[ "${1:-}" != "--" ]]; then
        echo "ERROR: expected -- before chromosome array" >&2
        return 1
    fi
    shift

    local chroms=()
    while [[ $# -gt 0 && "$1" != "--" ]]; do
        chroms+=( "$1" )
        shift
    done
    if [[ $# -eq 0 ]]; then
        echo "ERROR: expected -- before population array" >&2
        return 1
    fi
    shift
    local pops=( "$@" )
    local num_chroms="${#chroms[@]}"
    local num_pops="${#pops[@]}"
    local max_task=$((num_chroms * num_pops))
    local chr_index
    local pop_index

    if (( num_chroms < 1 )); then
        echo "ERROR: at least one chromosome is required" >&2
        return 1
    fi
    if (( num_pops < 1 )); then
        echo "ERROR: at least one population is required" >&2
        return 1
    fi
    if (( task_id < 1 || task_id > max_task )); then
        echo \
            "ERROR: task ID ${task_id} outside valid range 1..${max_task}" \
            >&2
        return 1
    fi

    # each chromosome owns a contiguous, population-minor task block.
    chr_index=$(((task_id - 1) / num_pops))
    pop_index=$(((task_id - 1) % num_pops))
    printf "%s %s\n" \
        "${chroms[chr_index]}" \
        "${pops[pop_index]}"
}
