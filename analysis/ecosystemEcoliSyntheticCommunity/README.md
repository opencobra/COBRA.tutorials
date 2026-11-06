# Abundance-growth space of a synthetic community of two _E. coli_ auxotrophs

| File | Description |
|---|---|
| [tutorial_ecosystemEcoliSyntheticCommunity.mlx](tutorial_ecosystemEcoliSyntheticCommunity.mlx) | MATLAB Live Script, saved with its outputs and figures |
| [tutorial_ecosystemEcoliSyntheticCommunity.m](tutorial_ecosystemEcoliSyntheticCommunity.m) | MATLAB script (generated from the Live Script) |

Reproduction of the notebook `sup_coli_synthetic_com.ipynb` of Jiménez et al. [1]
(https://github.com/mathomics/ecosystem): a community of a lysine auxotroph and a leucine auxotroph
of _E. coli_ (model _i_AF1260) that cross-feed the amino acid that the partner cannot synthesize.

Steps: construction of the auxotrophic strains and of the community model, selection of the
non-blocked reactions, abundance-growth grid and feasibility, qualitative flux variability analysis
(~3,000 reactions at 23 points), clustering of the grid points, reactions that characterize the
clusters and quantitative flux coupling between the lysine exchanges of both strains.

Requirements: MATLAB (no additional MATLAB toolbox is required), the COBRA Toolbox (functions in
`src/analysis/ecosystem`; the model `iAF1260.mat` is distributed with the toolbox) and an LP solver
(Gurobi or IBM CPLEX strongly recommended). Optionally, a parallel pool (Parallel Computing Toolbox)
speeds up the flux variability analyses.

Runtime (Gurobi 9.5 with a parallel pool of 3 workers, MATLAB R2021a): ~1.5 h, of which most
corresponds to the qualitative flux variability analysis.

[1] Jiménez, N.E., Acuña, V., Cortés, M.P., Eveillard, D., Maass, A. Unveiling abundance-dependent
metabolic phenotypes of microbial communities. *mSystems* **8**, e00492-23 (2023).
