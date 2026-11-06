# Abundance-growth space of a bioleaching community: analysis of multiple scenarios

| File | Description |
|---|---|
| [tutorial_ecosystemBioleachingScenarios.mlx](tutorial_ecosystemBioleachingScenarios.mlx) | MATLAB Live Script, saved with its outputs and figures |
| [tutorial_ecosystemBioleachingScenarios.m](tutorial_ecosystemBioleachingScenarios.m) | MATLAB script (generated from the Live Script) |

Reproduction of the notebook `sup_mat_multiple_bioleaching.ipynb` of Jiménez et al. [1]
(https://github.com/mathomics/ecosystem): a bioleaching community of _Acidithiobacillus
ferrooxidans_ Wenelen (_i_ML510) and _Sulfobacillus thermosulfidooxidans_ Cutipay (_i_SM517) simulated
under several scenarios of substrate supply, energy level, biological constraints and carbon
release as exopolysaccharides.

Steps: load and filter the table of simulated scenarios, load the results of each scenario (Python
pickle, NumPy and cobrapy JSON files), clusters of one and of all the scenarios, flux values of a
reaction on the grid, overlap of qualitative states among scenarios, key reaction candidates and
cluster descriptors.

Data: folder `result_files` and file `simulations_pareto4.csv` of https://github.com/mathomics/ecosystem.
Set the variable `ecosystemRepoDir` to a local copy of the repository before running the tutorial;
otherwise, the version of the repository used in the tutorial is downloaded (~35 MB) to a
temporary folder.

Requirements: MATLAB (no additional MATLAB toolbox is required), the COBRA Toolbox (functions in
`src/analysis/ecosystem` and `src/base/io/utilities`) and an LP solver (only to recompute the
feasibility of one scenario and the maximum flux of a reaction in Step 5).

Runtime: ~2 min (Gurobi 9.5, MATLAB R2021a), plus the download of the data.

[1] Jiménez, N.E., Acuña, V., Cortés, M.P., Eveillard, D., Maass, A. Unveiling abundance-dependent
metabolic phenotypes of microbial communities. *mSystems* **8**, e00492-23 (2023).
