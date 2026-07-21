###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_local_ancestry.py
###############################################################################
# write local ancestry intervals for downstream statistics.


##### set up ##################################################################
from .build_sample_ids import build_sample_ids


POP_ID_TO_NAME = {
    0: "AFR",
    1: "EUR",
    2: "ADX",
}


##### main function ###########################################################
'''
write tspop local ancestry intervals with chromosome, sample, haplotype, and
ancestral-population labels. The output omits a header for downstream tools.
'''
def write_local_ancestry(
    ancestry_table,
    sample_node_rows,
    out_tsv,
    chr_used,
    pop,
):
    sample_to_ind = {node: ind_id for node, ind_id, _ in sample_node_rows}
    sample_to_hap = {node: hap for node, _, hap in sample_node_rows}

    ancestry_table = ancestry_table.copy()
    pop_start_ind = min(sample_to_ind[int(node)] for node in ancestry_table[
        "sample"
    ])
    ancestry_table["sample_id"] = ancestry_table["sample"].map(
        lambda node: build_sample_ids(
            pop,
            sample_to_ind[int(node)],
            pop_start_ind,
        )[0]
    )
    ancestry_table["vcf_sample_id"] = ancestry_table["sample"].map(
        lambda node: build_sample_ids(
            pop,
            sample_to_ind[int(node)],
            pop_start_ind,
        )[1]
    )
    ancestry_table["hap"] = ancestry_table["sample"].map(
        lambda node: sample_to_hap[int(node)]
    )
    ancestry_table["population_name"] = ancestry_table["population"].map(
        lambda pop: POP_ID_TO_NAME.get(int(pop), f"pop_{int(pop)}")
    )
    ancestry_table.insert(0, "chrom", f"chr{chr_used}")
    ancestry_table = ancestry_table[
        [
            "chrom",
            "left",
            "right",
            "sample_id",
            "vcf_sample_id",
            "hap",
            "population_name",
            "population",
            "ancestor",
        ]
    ]
    ancestry_table.to_csv(out_tsv, sep="\t", header=False, index=False)
