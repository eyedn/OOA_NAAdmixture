###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_multik_ancestry_rows.py
###############################################################################

# overview: build neutral multi-K ancestry rows in authoritative sample order.


##### set up ##################################################################
from .normalize_multik_admixture_rows import (
    normalize_multik_admixture_rows,
)


##### main function ###########################################################
'''
normalize raw component matrices and return K-major, FAM-minor ancestry rows.
'''
def build_multik_ancestry_rows(
    samples,
    sample_pops,
    values_by_k,
    chrom=None,
    rep=0,
    sample_ids=None,
):
    output_rows = []
    for k in sorted(values_by_k):
        values = values_by_k[k]
        if len(values) != len(samples):
            raise ValueError(f"K={k} component/FAM row count mismatch")
        raw_rows = [
            {
                "sample": sample,
                "pop": sample_pops[sample],
                "q_values": q_values,
            }
            for sample, q_values in zip(samples, values)
        ]
        for normalized in normalize_multik_admixture_rows(raw_rows, k):
            row = {"rep": rep}
            if chrom is not None:
                row["chrom"] = chrom
            row.update(
                {
                    "pop": normalized["pop"],
                    "sample_id": (
                        sample_ids[normalized["sample"]]
                        if sample_ids is not None
                        else "NA"
                    ),
                    "vcf_sample_id": normalized["sample"],
                    "k": k,
                    "component_1_q": normalized["component_1_q"],
                    "component_2_q": normalized["component_2_q"],
                    "component_3_q": normalized["component_3_q"],
                    "component_4_q": normalized["component_4_q"],
                    "component_5_q": normalized["component_5_q"],
                    "span": "NA",
                }
            )
            output_rows.append(row)
    return output_rows
