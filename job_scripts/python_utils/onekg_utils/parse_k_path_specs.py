###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           parse_k_path_specs.py
###############################################################################

# overview: validate repeated K=PATH inference-output specifications.


##### main function ###########################################################
'''
return requested inference paths keyed by K after validating each specification.
'''
def parse_k_path_specs(specs, tool_name):
    paths = {}
    for spec in specs:
        k_text, separator, path = spec.partition("=")
        if not separator or not k_text.isdigit() or not path:
            raise ValueError(
                f"{tool_name} paths must use the format K=PATH"
            )
        k = int(k_text)
        if k < 2 or k > 5:
            raise ValueError(f"{tool_name} K must be between 2 and 5")
        if k in paths:
            raise ValueError(f"Duplicate {tool_name} K={k}")
        paths[k] = path
    if not paths:
        raise ValueError(f"At least one {tool_name} K=PATH is required")
    return paths
