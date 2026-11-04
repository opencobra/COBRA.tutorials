# Abundance-growth space of a microbial community: toy example

| File | Description |
|---|---|
| [tutorial_ecosystemToyExample.mlx](tutorial_ecosystemToyExample.mlx) | MATLAB Live Script, saved with its outputs and figures |
| [tutorial_ecosystemToyExample.m](tutorial_ecosystemToyExample.m) | MATLAB script (generated from the Live Script) |

Reproduction of the notebook `sup_toy_example.ipynb` of Jiménez et al. [1]
(https://github.com/mathomics/ecosystem): a community of two organisms that compete for substrate A,
while organism 2 produces B, which can be used by organism 1.

Steps: build the member models and the community model, set the simulation conditions, select the
reactions to analyze, build the abundance-growth grid, check the feasibility of the grid points,
qualitative flux variability analysis, clustering of the grid points, representative qualitative
states of each cluster and quantitative flux coupling analysis.

Requirements: MATLAB (no additional MATLAB toolbox is required), the COBRA Toolbox (functions in
`src/analysis/ecosystem`) and an LP solver.

Runtime: ~3 min (Gurobi 9.5, MATLAB R2021a).

[1] Jiménez, N.E., Acuña, V., Cortés, M.P., Eveillard, D., Maass, A. Unveiling abundance-dependent
metabolic phenotypes of microbial communities. *mSystems* **8**, e00492-23 (2023).
