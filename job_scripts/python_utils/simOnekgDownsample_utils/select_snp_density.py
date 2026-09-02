###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           select_snp_density.py
###############################################################################

# overview: select simulation SNPs to match empirical counts by fixed window.


##### set up ##################################################################
from random import Random


##### main ####################################################################
'''
return selected positions and union-window density contract rows.
'''
def select_snp_density(
    empirical_rows,
    simulation_rows,
    simulation_positions,
    window_size_bp,
    seed,
    rep
):
    if window_size_bp <= 0:
        raise ValueError("Window size must be positive")

    empirical_by_window = {}
    simulation_by_window = {}
    chromosomes = set()
    for rows, output in (
        (empirical_rows, empirical_by_window),
        (simulation_rows, simulation_by_window)
    ):
        for row in rows:
            chrom = str(row["chrom"])
            chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
            bin_start = int(row["bin_start"])
            count = int(row["snp_count"])
            key = (chrom, bin_start)
            if key in output:
                raise ValueError(f"Duplicate density window: {key}")
            if bin_start < 0 or count < 0:
                raise ValueError(
                    "Density coordinates and counts must be nonnegative"
                )
            output[key] = count
            chromosomes.add(chrom)

    positions_by_window = {}
    seen_positions = set()
    source_chrom_by_canonical = {}
    for raw_chrom, raw_position in simulation_positions:
        source_chrom = str(raw_chrom)
        chrom = (
            source_chrom[3:]
            if source_chrom.lower().startswith("chr")
            else source_chrom
        )
        position = int(raw_position)
        key = (chrom, position)
        if key in seen_positions:
            raise ValueError(
                f"Duplicate simulation position: {source_chrom}:{position}"
            )
        if position <= 0:
            raise ValueError(
                "Simulation positions must be one-based and positive"
            )
        seen_positions.add(key)
        source_chrom_by_canonical[chrom] = source_chrom
        # Match vcftools --SNPdensity exactly: idx = POS / bin_size.
        bin_start = (position // window_size_bp) * window_size_bp
        positions_by_window.setdefault((chrom, bin_start), []).append(position)
        chromosomes.add(chrom)

    if len(chromosomes) != 1:
        raise ValueError(
            "Density inputs and positions must describe one chromosome"
        )
    for key, positions in positions_by_window.items():
        expected_count = simulation_by_window.get(key, 0)
        if len(positions) != expected_count:
            raise ValueError(
                f"Simulation density count disagrees with positions for {key}: "
                f"density count={expected_count}, positions={len(positions)}"
            )
    for key, expected_count in simulation_by_window.items():
        if len(positions_by_window.get(key, [])) != expected_count:
            raise ValueError(
                f"Simulation density count disagrees with positions for {key}: "
                f"density count={expected_count}, "
                f"positions={len(positions_by_window.get(key, []))}"
            )

    rng = Random(seed)
    selected = []
    contracted_rows = []
    all_windows = sorted(set(empirical_by_window) | set(simulation_by_window))
    for chrom, bin_start in all_windows:
        empirical_count = empirical_by_window.get((chrom, bin_start), 0)
        simulation_count = simulation_by_window.get((chrom, bin_start), 0)
        selected_count = min(empirical_count, simulation_count)
        available = positions_by_window.get((chrom, bin_start), [])
        selected_positions = (
            available
            if selected_count == len(available)
            else rng.sample(available, selected_count)
        )
        source_chrom = source_chrom_by_canonical.get(chrom, chrom)
        selected.extend(
            (source_chrom, position) for position in selected_positions
        )
        deficit = max(empirical_count - simulation_count, 0)
        rate_divisor = window_size_bp / 1000
        contracted_rows.append(
            {
                "rep": rep,
                "chrom": chrom,
                "bin_start": bin_start,
                "bin_end": bin_start + window_size_bp,
                "window_size_bp": window_size_bp,
                "empirical_snp_count": empirical_count,
                "simulation_snp_count": simulation_count,
                "selected_snp_count": selected_count,
                "empirical_variants_per_kb": empirical_count / rate_divisor,
                "simulation_variants_per_kb": simulation_count / rate_divisor,
                "selected_variants_per_kb": selected_count / rate_divisor,
                "deficit_snp_count": deficit,
                "simulation_below_target": deficit > 0
            }
        )

    selected.sort(
        key=lambda value: (
            value[0][3:] if value[0].lower().startswith("chr") else value[0],
            value[1]
        )
    )
    return selected, contracted_rows
