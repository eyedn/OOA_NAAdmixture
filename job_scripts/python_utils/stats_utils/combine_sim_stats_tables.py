###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_sim_stats_tables.py
###############################################################################


from pathlib import Path
import pandas as pd


# combine per-replicate statistics tables after the stats array succeeds
def combine_sim_stats_tables(args, table_names):
    stats_path = Path(args.stats_dir)
    for table_name in table_names:
        per_rep_paths = [
            stats_path / f"{table_name}.rep_{rep}.tsv"
            for rep in range(1, args.num_reps + 1)
        ]
        missing_paths = [path for path in per_rep_paths if not path.exists()]
        if missing_paths:
            missing = ", ".join(str(path) for path in missing_paths)
            raise FileNotFoundError(
                f"Missing per-replicate files for {table_name}: {missing}"
            )

        combined = pd.concat(
            [pd.read_csv(path, sep="\t") for path in per_rep_paths],
            ignore_index=True,
        )
        combined.to_csv(
            stats_path / f"{table_name}.tsv",
            sep="\t",
            index=False,
        )
        combined.to_parquet(
            stats_path / f"{table_name}.parquet",
            index=False,
        )
