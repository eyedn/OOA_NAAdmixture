###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_local_ancestry.py
###############################################################################

# pattern: Imperative Shell


from .build_adx_sample_ids import build_adx_sample_ids


POP_ID_TO_NAME = {
    0: "AFR",
    1: "EUR",
    2: "ADX",
}


# output tsv of tspop local ancestry
def write_local_ancestry(ancestry_table, sample_node_rows, out_tsv, chr_used):
    sample_to_ind = {node: ind_id for node, ind_id, _ in sample_node_rows}
    sample_to_hap = {node: hap for node, _, hap in sample_node_rows}

    ancestry_table = ancestry_table.copy()
    adx_start_ind = min(sample_to_ind[int(node)] for node in ancestry_table[
        "sample"
    ])
    ancestry_table["sample_id"] = ancestry_table["sample"].map(
        lambda node: build_adx_sample_ids(
            sample_to_ind[int(node)],
            adx_start_ind,
        )[0]
    )
    ancestry_table["vcf_sample_id"] = ancestry_table["sample"].map(
        lambda node: build_adx_sample_ids(
            sample_to_ind[int(node)],
            adx_start_ind,
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
