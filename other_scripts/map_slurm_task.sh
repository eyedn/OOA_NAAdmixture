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


# map a slurm array task id to a replicate and chromosome
map_slurm_task_to_rep_chr() {
    local task_id="$1"
    local num_reps="$2"
    shift 2
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

    chr_index=$(((task_id - 1) / num_reps))
    rep=$((((task_id - 1) % num_reps) + 1))
    chr="${chroms[chr_index]}"

    printf "%s %s\n" "${rep}" "${chr}"
}
