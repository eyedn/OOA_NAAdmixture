###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           project_complete_sfs.py
###############################################################################

# overview: project raw sfs counts to a downsampled size.


##### set up ##################################################################
from functools import lru_cache
from math import comb


##### internal functions ######################################################
'''
internal: return the weights that dictate the downsampling of sfs bins.
'''
@lru_cache(maxsize=None)
def _build_projection_weights(source_size, target_size):
    # return dense expected hypergeometric projection weights.
    denominator = comb(source_size, target_size)
    weights = []
    for target_count in range(target_size + 1):
        row = []
        for source_count in range(source_size + 1):
            if (
                target_count > source_count
                or target_size - target_count > source_size - source_count
            ):
                row.append(0.0)
                continue
            # expected resampling projection; see Mah et al. (2025),
            # Gutenkunst et al. (2009), and the dadi projection documentation.
            numerator = comb(source_count, target_count) * comb(
                source_size - source_count,
                target_size - target_count
            )
            row.append(numerator / denominator)
        weights.append(tuple(row))
    return tuple(weights)


##### main function ###########################################################
'''
project a complete unfolded spectrum to ``target_size`` copies.
'''
def project_complete_sfs(spectrum, target_size):
    values = tuple(float(value) for value in spectrum)
    source_size = len(values) - 1
    weights = _build_projection_weights(source_size, target_size)
    return tuple(
        sum(weight * value for weight, value in zip(row, values))
        for row in weights
    )
