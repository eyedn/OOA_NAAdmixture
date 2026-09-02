###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           parse_sim_positions.py
###############################################################################

# overview: normalize CHROM/POS rows queried from a simulation VCF.


##### main ####################################################################
'''
parse two-column CHROM/POS text rows.
'''
def parse_sim_positions(lines):
    positions = []
    for line in lines:
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 2:
            raise ValueError("Position rows must contain exactly CHROM and POS")
        positions.append((fields[0], int(fields[1])))
    if not positions:
        raise ValueError("Simulation position table is empty")
    return positions
