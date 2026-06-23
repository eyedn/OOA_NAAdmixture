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
    "build_demography": (".build_demography", "build_demography"),
    "build_sample_metadata_rows": (
        ".build_sample_metadata_rows",
        "build_sample_metadata_rows",
    ),
    "run_simulation": (".run_simulation", "run_simulation"),
}


def __getattr__(name):
    from importlib import import_module

    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    return getattr(import_module(module_name, __name__), attr_name)


__all__ = [
    "build_demography",
    "build_sample_metadata_rows",
    "run_simulation",
]
