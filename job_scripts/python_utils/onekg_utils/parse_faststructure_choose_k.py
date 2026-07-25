###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           parse_faststructure_choose_k.py
###############################################################################

# overview: validate and normalize a fastStructure chooseK report.


##### set up ##################################################################
import re


REPORT_PATTERNS = {
    "max_marginal_likelihood_k": re.compile(
        r"^\s*Model complexity that maximizes marginal likelihood"
        r"\s*=\s*(\S+)\s*$",
        re.MULTILINE,
    ),
    "model_components_k": re.compile(
        r"^\s*Model components used to explain structure in data"
        r"\s*=\s*(\S+)\s*$",
        re.MULTILINE,
    ),
}


##### main function ###########################################################
'''
parse the two chooseK recommendations and attach inference provenance.
'''
def parse_faststructure_choose_k(
    report,
    prior,
    seed,
    requested_ks,
    chrom=None,
):
    ks = sorted(set(requested_ks))
    if not ks or len(ks) != len(requested_ks):
        raise ValueError("fastStructure K values must be nonempty and unique")
    parsed = {}
    for name, pattern in REPORT_PATTERNS.items():
        matches = pattern.findall(report)
        if len(matches) != 1:
            raise ValueError(
                f"chooseK report must contain exactly one {name} value"
            )
        try:
            parsed[name] = int(matches[0])
        except ValueError as error:
            raise ValueError(
                f"chooseK {name} must be an integer"
            ) from error
        if parsed[name] < ks[0] or parsed[name] > ks[-1]:
            raise ValueError(
                f"chooseK {name} is outside requested K range"
            )
    row = {"rep": 0}
    if chrom is not None:
        row["chrom"] = chrom
    row.update(
        {
            "prior": prior,
            "seed": seed,
            "k_min": ks[0],
            "k_max": ks[-1],
            **parsed,
        }
    )
    return row
