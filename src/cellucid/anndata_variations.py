import gc
import pandas as pd
import itertools
from pathlib import Path
from tqdm import tqdm

try:
    import scanpy as sc
except ImportError as exc:  # pragma: no cover - optional dependency
    raise ImportError(
        "scanpy is required for cellucid.anndata_variations; install with "
        "`pip install cellucid[analysis]` or add scanpy to your environment."
    ) from exc

def _find_repo_root() -> Path:
    """Locate the python project root (contains pyproject.toml and data/)."""
    for candidate in Path(__file__).resolve().parents:
        if (candidate / "pyproject.toml").exists() and (candidate / "data").exists():
            return candidate
    return Path.cwd()


REPO_ROOT = _find_repo_root()
RAW_DATA_DIR = REPO_ROOT / "data" / "raw"
EXPERIMENT_DIR = REPO_ROOT / "data" / "experiments"

INPUT_H5AD = RAW_DATA_DIR / "adata_latent_scvi-adata_unified_20250925_001_integration-supergpu02.scidom.de-bb886d36409e78a3.h5ad"
METADATA_PICKLE = RAW_DATA_DIR / "adata_obs_extracted.pickle"
OUTPUT_DIR = EXPERIMENT_DIR / "umap_parameter_sweep"

MIN_DISTS = [0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
N_NEIGHBORS_LIST = [5, 10, 15, 30, 50, 90, 120]
RANDOM_SEED = 42


def _print_setup():
    """Log resolved paths so it's obvious what is being used."""
    print("==== AnnData UMAP sweep setup ====")
    print(f"Current working dir : {Path.cwd()}")
    print(f"Repo root           : {REPO_ROOT}")
    print(f"Input h5ad          : {INPUT_H5AD} (exists: {INPUT_H5AD.exists()})")
    print(f"Metadata pickle     : {METADATA_PICKLE} (exists: {METADATA_PICKLE.exists()})")
    print(f"Output directory    : {OUTPUT_DIR}")
    print("==================================")


def validate_inputs() -> None:
    """Ensure required files and folders exist before heavy work."""
    missing = [path for path in (INPUT_H5AD, METADATA_PICKLE) if not path.exists()]
    if missing:
        missing_list = "\n - ".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing required input files:\n - {missing_list}")

    EXPERIMENT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)



def load_and_preprocess_data():
    validate_inputs()
    print("Loading and preparing data (this happens only once)...")
    adata = sc.read_h5ad(INPUT_H5AD)
    adata_complete = pd.read_pickle(METADATA_PICKLE)

    if "qc_fail" not in adata_complete.columns:
        raise KeyError("Metadata pickle missing required column 'qc_fail'.")

    missing_cells = set(adata.obs_names) - set(adata_complete.index)
    if missing_cells:
        raise ValueError(f"Metadata missing {len(missing_cells)} cells (e.g., {next(iter(missing_cells))}).")

    adata_complete = adata_complete.loc[adata.obs_names]
    qc_mask = ~adata_complete["qc_fail"].fillna(False)

    # Apply QC filtering to both the matrix and obs
    adata = adata[qc_mask].copy()
    adata_complete = adata_complete.loc[qc_mask]
    adata.obs = adata_complete

    # Drop heavy derived fields to reduce per-run memory
    for container in (adata.obsm, adata.obsp):
        for key in list(container.keys()):
            del container[key]
    adata.uns.clear()

    print(f"Data prepared. Final: {adata}")
    print(f"Cells kept after QC : {adata.n_obs}")
    print(f"Features            : {adata.n_vars}")
    return adata


def run_parameter_sweep(_adata):
    """
    Iterates through parameters, calculating UMAPs and saving results.
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    param_combinations = list(itertools.product(N_NEIGHBORS_LIST, MIN_DISTS))
    print(f"Starting sweep over {len(param_combinations)} combinations...")
    print(f"Writing outputs to: {OUTPUT_DIR}")

    existing_before = 0
    for k, min_dist in param_combinations:
        dist_str = f"{min_dist:.2f}".replace('.', '_')
        filename = OUTPUT_DIR / f"adata|k__{k}|min_dist_{dist_str}.h5ad"
        if filename.exists():
            existing_before += 1
    to_compute = len(param_combinations) - existing_before
    print(f"Found {existing_before} outputs already present; {to_compute} will be computed.")

    if to_compute == 0:
        print("All combinations already exist. Nothing to do.")
        return

    skipped = 0
    created = 0

    for k, min_dist in tqdm(param_combinations, desc="Processing"):
        dist_str = f"{min_dist:.2f}".replace('.', '_')
        filename = OUTPUT_DIR / f"adata|k__{k}|min_dist_{dist_str}.h5ad"
        if filename.exists():
            skipped += 1
            tqdm.write(f"Exists, skipping: {filename.name}")
            continue

        adata = _adata.copy()
        try:
            sc.pp.neighbors(adata, n_neighbors=k, random_state=RANDOM_SEED)
            sc.tl.umap(
                adata, 
                n_components=3, 
                min_dist=min_dist, 
                random_state=RANDOM_SEED,
            )
            adata.write_h5ad(filename)
            created += 1
        except Exception as e:
            tqdm.write(f"Calculation error for {filename.name}: {e}")
        gc.collect()

    print(f"Finished sweep. Created {created}, skipped {skipped}, total {len(param_combinations)}.")

if __name__ == "__main__":
    try:
        _print_setup()
        clean_adata = load_and_preprocess_data()
        run_parameter_sweep(clean_adata)
        print("Processing complete.")
    except Exception as e:
        print(f"An error occurred: {e}")
