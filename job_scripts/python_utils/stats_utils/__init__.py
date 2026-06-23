###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           __init__.py
###############################################################################

# pattern: Imperative Shell


_EXPORTS = {
    "build_1d_sfs_rows": (".build_1d_sfs_rows", "build_1d_sfs_rows"),
    "build_2d_sfs_rows": (".build_2d_sfs_rows", "build_2d_sfs_rows"),
    "build_ld_decay_rows": (".build_ld_decay_rows", "build_ld_decay_rows"),
    "build_pi_theta_rows": (".build_pi_theta_rows", "build_pi_theta_rows"),
    "calc_stats": (".calc_stats", "calc_stats"),
    "combine_sim_stats_tables": (
        ".combine_sim_stats_tables",
        "combine_sim_stats_tables",
    ),
    "parse_admixture_q": (".parse_admixture_q", "parse_admixture_q"),
    "parse_king_file": (".parse_king_file", "parse_king_file"),
    "read_fam_order": (".read_fam_order", "read_fam_order"),
    "read_sample_metadata": (".read_sample_metadata", "read_sample_metadata"),
    "sample_nodes_for_pop": (".sample_nodes_for_pop", "sample_nodes_for_pop"),
    "SimpleTable": (".simple_table", "SimpleTable"),
    "subset_pairs": (".subset_pairs", "subset_pairs"),
    "write_admixture_pop": (".write_admixture_pop", "write_admixture_pop"),
}


def __getattr__(name):
    from importlib import import_module

    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    return getattr(import_module(module_name, __name__), attr_name)


__all__ = [
    "build_1d_sfs_rows",
    "build_2d_sfs_rows",
    "build_ld_decay_rows",
    "build_pi_theta_rows",
    "calc_stats",
    "combine_sim_stats_tables",
    "parse_admixture_q",
    "parse_king_file",
    "read_fam_order",
    "read_sample_metadata",
    "sample_nodes_for_pop",
    "SimpleTable",
    "subset_pairs",
    "write_admixture_pop",
]
