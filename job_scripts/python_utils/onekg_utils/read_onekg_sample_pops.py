###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           read_onekg_sample_pops.py
###############################################################################

# overview: map unrelated 1000 Genomes samples to requested populations.


##### main function ###########################################################
'''
read requested sample populations and validate them against a PLINK FAM file.
Returns a sample-to-population mapping limited to requested populations.
'''
def read_onekg_sample_pops(unrels_path, fam_path, pops):
    with open(unrels_path, "r", encoding="utf-8") as in_file:
        lines = [line.split() for line in in_file if line.strip()]
    if not lines:
        raise ValueError("ONEKG_UNRELS_FILE is empty")
    headers = [header.lower() for header in lines[0]]
    iid_names = ("iid", "sample", "sample_id", "individual", "id")
    pop_names = ("pop", "population", "population_code")
    iid_index = next(
        (headers.index(name) for name in iid_names if name in headers),
        None,
    )
    pop_index = next(
        (headers.index(name) for name in pop_names if name in headers),
        None,
    )
    if iid_index is None or pop_index is None:
        raise ValueError(
            "Unrelated-sample file needs IID and population columns"
        )
    sample_pops = {}
    for fields in lines[1:]:
        if max(iid_index, pop_index) >= len(fields):
            raise ValueError("Malformed row in unrelated-sample file")
        sample = fields[iid_index]
        pop = fields[pop_index]
        if pop in pops:
            if sample in sample_pops and sample_pops[sample] != pop:
                raise ValueError(f"Conflicting populations for IID {sample}")
            sample_pops[sample] = pop
    if not sample_pops:
        raise ValueError(
            "No requested populations found in unrelated-sample file"
        )

    with open(fam_path, "r", encoding="utf-8") as in_file:
        fam_samples = {
            fields[1]
            for line in in_file
            if len(fields := line.split()) >= 2
        }
    missing = sorted(set(sample_pops) - fam_samples)
    if missing:
        raise ValueError(
            "Requested IIDs missing from configured FAM: " + ", ".join(missing)
        )
    return sample_pops
