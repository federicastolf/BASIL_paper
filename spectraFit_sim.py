
!pip install scanpy anndata
!pip install --no-deps git+https://github.com/dpeerlab/spectra.git
!pip install pyvis

import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from Spectra import Spectra

!unzip -q sim_data.zip

SIM_DIR = "sim_data"
OUT_DIR = "spectra_results"
os.makedirs(OUT_DIR, exist_ok=True)

sim_dirs = sorted([
    d for d in os.listdir(SIM_DIR)
    if d.startswith("sim_") and os.path.isdir(os.path.join(SIM_DIR, d))
])

import Spectra.Spectra as spectra_module
import numpy as np
import pandas as pd

def return_markers_fixed(factor_matrix, id2word, n_top_vals=100):
    idx_matrix = np.argsort(factor_matrix, axis=1)[:, ::-1][:, :n_top_vals]
    df = pd.DataFrame(np.empty(idx_matrix.shape, dtype=object))
    for i in range(idx_matrix.shape[0]):
        for j in range(idx_matrix.shape[1]):
            df.iloc[i, j] = id2word[idx_matrix[i, j]]
    return df.values

spectra_module.return_markers = return_markers_fixed

for sim_id in sim_dirs:
    sim_path = os.path.join(SIM_DIR, sim_id)
    out_path = os.path.join(OUT_DIR, sim_id)
    os.makedirs(out_path, exist_ok=True)

    print(f"[{sim_id}] loading...")

    # --- load data ---
    Y = np.loadtxt(os.path.join(sim_path, "Y.csv"), delimiter=",")   # n x p
    C = np.loadtxt(os.path.join(sim_path, "C.csv"), delimiter=",")   # p x q
    with open(os.path.join(sim_path, "gene_set_names.txt")) as f:
        gs_names = [line.strip() for line in f if line.strip()]

    n, p = Y.shape
    q = C.shape[1]
    assert C.shape[0] == p
    assert len(gs_names) == q
    gene_names = [f"g{j}" for j in range(p)]

    # build AnnData
    adata = ad.AnnData(
        X=Y.astype(np.float32),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=gene_names)
    )

    # build gene_set_dictionary
    gene_set_dictionary = {}
    for j, gs in enumerate(gs_names):
        member_idx = np.where(C[:, j] > 0)[0]
        if len(member_idx) > 0:
            gene_set_dictionary[gs] = [gene_names[i] for i in member_idx]

    print(f"[{sim_id}] fitting Spectra: n={n}, p={p}, q={len(gene_set_dictionary)}")

    # --- fit ---
    t0 = time.time()
    model = Spectra.est_spectra(
        adata=adata,
        gene_set_dictionary=gene_set_dictionary,
        L=10,
        use_cell_types=False,
        use_highly_variable=False,
        use_weights=True,
        lam=0.1,
        delta=0.001,
        kappa=None,
        rho=0.001,
        n_top_vals=50,
        label_factors=False,
        num_epochs=5000,
    )
    elapsed = time.time() - t0

    # extract factor and loadings
    Theta = np.asarray(adata.uns["SPECTRA_factors"])       # K x p
    A     = np.asarray(adata.obsm["SPECTRA_cell_scores"])  # n x K

    Lambda_hat = Theta.T   # p x K
    M_hat = A              # n x K

    # save results
    np.savetxt(os.path.join(out_path, "Lambda_hat.csv"), Lambda_hat, delimiter=",")
    np.savetxt(os.path.join(out_path, "M_hat.csv"),      M_hat,      delimiter=",")
    with open(os.path.join(out_path, "meta.txt"), "w") as f:
        f.write(f"time_seconds\t{elapsed}\n")
        f.write(f"K\t{Theta.shape[0]}\n")

    print(f"[{sim_id}] done in {elapsed:.1f}s, K={Theta.shape[0]}")

!zip -r spectra_results.zip spectra_results

# from google.colab import files
# files.download("spectra_results.zip")