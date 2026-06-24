###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           calc_stats.py
###############################################################################


from pathlib import Path
import pandas as pd
import tszip
from .build_2d_sfs_rows import build_2d_sfs_rows
from .build_ancestry_table import build_ancestry_table
from .build_ld_decay_rows import build_ld_decay_rows
from .build_pi_theta_rows import build_pi_theta_rows
from .build_1d_sfs_rows import build_1d_sfs_rows
from .parse_king_file import parse_king_file


# generate replicate summaries on pi, theta, sfs, ld, and kinship
def calc_stats(args):
    stats_dir = Path(args.stats_dir)
    stats_dir.mkdir(parents=True, exist_ok=True)

    # load compressed tree sequence
    tree_tsz_path = Path(args.tree_tsz_path)
    if not str(tree_tsz_path).endswith(".ts.tsz"):
        raise ValueError(
            f"Tree sequence input must end with .ts.tsz: {tree_tsz_path}"
        )
    if not tree_tsz_path.exists():
        raise FileNotFoundError(
            f"Missing compressed tree sequence {tree_tsz_path}"
        )
    ts = tszip.decompress(tree_tsz_path)

    # process tspop and "admixture --supervised" ancestry results
    q_path = Path(args.admixture_dir) / (
        f"{args.genetic_map}_{args.rep}_all.2.Q"
    )
    ancestry_table = build_ancestry_table(
        rep=args.rep,
        pops=args.pops,
        genetic_map=args.genetic_map,
        global_anc_dir=args.global_anc_dir,
        q_path=q_path,
        fam_path=args.admixture_fam_path,
    )
    ancestry_table.to_csv(
        stats_dir / f"ancestry.rep_{args.rep}.tsv",
        sep="\t",
        index=False
    )
    ancestry_table.to_parquet(
        stats_dir / f"ancestry.rep_{args.rep}.parquet",
        index=False
    )

    # process KING kinship tables
    king_tables = []
    for pop in args.pops:
        king_path = Path(args.king_dir) / (
            f"{args.genetic_map}_{args.rep}_{pop}.kin0"
        )
        king_table = parse_king_file(king_path, args.rep, pop)
        king_table.to_csv(
            Path(args.king_dir) / f"{args.genetic_map}_{args.rep}_{pop}.tsv",
            sep="\t",
            index=False
        )
        king_tables.append(king_table)
    king_combined = pd.concat(king_tables, ignore_index=True)
    king_combined.to_csv(
        stats_dir / f"kinship.rep_{args.rep}.tsv",
        sep="\t",
        index=False
    )
    king_combined.to_parquet(
        stats_dir / f"kinship.rep_{args.rep}.parquet",
        index=False
    )

    # derive and process pi and theta from the tree sequence
    pi_theta = pd.DataFrame(
        build_pi_theta_rows(
            ts,
            args.rep,
            args.pops,
            args.sample_size,
            args.mutation_rate,
        )
    )

    # derive and process 1D and 2D sfs from the tree sequence
    sfs = pd.DataFrame(
        build_1d_sfs_rows(ts, args.rep, args.pops, args.sample_size)
    )
    sfs_2d = pd.DataFrame(
        build_2d_sfs_rows(ts, args.rep, args.pops, args.sample_size)
    )

    # derive and process ld decay from the tree sequence
    ld_decay = pd.DataFrame(
        build_ld_decay_rows(ts, args.rep, args.pops, args.sample_size)
    )

    # combine summary tables into a singular table 
    tables = {
        "pi_theta_stats": pi_theta,
        "sfs": sfs,
        "sfs_2d": sfs_2d,
        "ld_decay": ld_decay,
    }
    for table_name, table in tables.items():
        table.to_csv(
            stats_dir / f"{table_name}.rep_{args.rep}.tsv",
            sep="\t",
            index=False,
        )
        table.to_parquet(
            stats_dir / f"{table_name}.rep_{args.rep}.parquet",
            index=False,
        )
