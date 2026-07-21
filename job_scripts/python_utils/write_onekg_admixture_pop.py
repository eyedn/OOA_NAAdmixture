###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_onekg_admixture_pop.py
###############################################################################

# overview: write 1000 Genomes supervised ADMIXTURE population labels.


##### set up ##################################################################
import argparse
from onekg_utils.read_onekg_sample_pops import read_onekg_sample_pops


##### internal functions #####################################################
'''
internal: return labels in final FAM order, retaining both references as
supervised and the admixed population as unsupervised.
'''
def _build_onekg_admixture_labels(
    fam_samples,
    sample_pops,
    afr_pop,
    eur_pop,
    admixed_pop,
):
    missing = sorted(set(fam_samples) - set(sample_pops))
    if missing:
        raise ValueError(
            "ADMIXTURE FAM samples missing population labels: "
            + ", ".join(missing)
        )
    fam_pops = [sample_pops[sample] for sample in fam_samples]
    missing_refs = [pop for pop in (afr_pop, eur_pop) if pop not in fam_pops]
    if missing_refs:
        raise ValueError(
            "ADMIXTURE FAM is missing reference labels: "
            + ", ".join(missing_refs)
        )
    labels = {afr_pop: afr_pop, eur_pop: eur_pop, admixed_pop: "-"}
    unexpected = sorted(set(fam_pops) - set(labels))
    if unexpected:
        raise ValueError(
            "ADMIXTURE FAM has unexpected populations: "
            + ", ".join(unexpected)
        )
    return [labels[pop] for pop in fam_pops]


##### arguments ###############################################################
'''
define empirical metadata and FAM paths, output path, and population labels for
the supervised ADMIXTURE population file.
'''
parser = argparse.ArgumentParser()
parser.add_argument("--unrels-path", required=True)
parser.add_argument("--source-fam-path", required=True)
parser.add_argument("--admixture-fam-path", required=True)
parser.add_argument("--pop-path", required=True)
parser.add_argument("--afr-pop", required=True)
parser.add_argument("--eur-pop", required=True)
parser.add_argument("--admixed-pop", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    pops = [args.afr_pop, args.eur_pop, args.admixed_pop]
    sample_pops = read_onekg_sample_pops(
        args.unrels_path,
        args.source_fam_path,
        pops,
    )
    with open(args.admixture_fam_path, "r", encoding="utf-8") as in_file:
        fam_samples = [line.split()[1] for line in in_file if line.strip()]
    labels = _build_onekg_admixture_labels(
        fam_samples,
        sample_pops,
        args.afr_pop,
        args.eur_pop,
        args.admixed_pop,
    )
    with open(args.pop_path, "w", encoding="utf-8") as out_file:
        out_file.writelines(f"{label}\n" for label in labels)
