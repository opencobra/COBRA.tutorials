%% Creation and simulation of personalized microbiota models through metagenomic data integration
%% Author: Federico Baldini, Molecular Systems Physiology Group, University of Luxembourg.
% Almut Heinken, 12/2020: Rewrote the tutorial to streamline it and add new features.
%% INTRODUCTION
% This tutorial shows the steps that MgPipe automatically performs to create 
% and simulate personalized microbiota models through metagenomic data integration. 
% Please note that this tutorial uses as an example a small dataset (4 columns 
% and 30 rows) with the purpose of demonstrating the functionalities of the pipeline. 
% We recommend using high-performance computing clusters when assembling and simulating 
% from bigger datasets. 
% 
% The pipeline is divided into 3 parts:
% 
% # *[PART 1]* Analysis of individuals' specific microbes abundances is computed. 
% Individuals' metabolic diversity in relation to microbiota size and disease 
% presence, as well as, classical multidimensional scaling (PCoA) on individuals' 
% reaction repertoire are examples.
% # *[PART 2]*: 1 Constructing a global metabolic model (setup) containing all 
% the microbes listed in the study. 2 Building individuals' specific models integrating 
% abundance data retrieved from metagenomics. For each organism, reactions are 
% coupled to their objective function. 
% # *[PART 3]* A specific range of growth is imposed for each microbiota model 
% and Simulations under specific diet regimes are carried. Set of standard analysis 
% to apply to the personalized models. PCA of computed MNPCs of individuals as 
% for example.
%% USAGE
% Normally, once provided all the input variables in the driver (StartMgPipe), 
% the only action required is to run the driver itself. However, for this tutorial, 
% we will disable the autorun functionality and compute each section manually. 
%% DRIVER
% This file has to be modified by the user to launch the pipeline and to define 
% inputs and outputs files and locations. 
%% Initialize the COBRA Toolbox
%%
initCobraToolbox(false)
%% Prepare input data and models
% We first set the paths to input and output files change directory to where 
% the tutorial is located
%%
tutorialPath = fileparts(which('tutorial_mgPipe'));
cd(tutorialPath);
%% 
% We will use the AGORA resource (Magnusdottir et al., Nat Biotechnol. 2017 
% Jan;35(1):81-89) in this tutorial. AGORA version 1.03 is available at https://github.com/VirtualMetabolicHuman/AGORA. 
% Download AGORA and place the models into a folder.
%%
system('curl -LJO https://github.com/VirtualMetabolicHuman/AGORA/archive/master.zip')
unzip('AGORA-master')
modPath = [pwd filesep 'AGORA-master' filesep 'CurrentVersion' filesep 'AGORA_1_03' filesep' 'AGORA_1_03_mat'];

% path where to save results
mkdir('Results');
resPath = [tutorialPath filesep 'Results'];
%% 
% path to and name of the file with dietary information. Here, we will use 
% an "Average European" diet that is located in the DietImplementation folder.
%%
dietFilePath='AverageEuropeanDiet';
%% 
% Then we set the path and the name of the file from which to load the abundances. 
% For this tutorial, to reduce the time of computations, we will use a reduced 
% version of the example file (normCoverageReduced.csv) provided in the folder 
% Resources: only 4 individuals and 30 strains will be considered. Plese, note 
% that abundances are normalized to a total sum of one. 
% 
% 
%%
abunFilePath='normCoverageReduced.csv';
%% 
% Next inputs will define:
% 
% # name of the objective function of organisms
% # format to use to save images
% # number of cores to use for the pipeline execution 
% # if to enable automatic detection and correction of possible bugs
% # if to enable compatibility mode 
% # if stratification criteria are available
% # if to simulate also a rich diet
% # if to use an external solver and save models with diet
% # the type of FVA function to use to solve 
% 
% The following setting should work for almost any system, but please check 
% carefully to be sure these options are valid for you. A more detailed description 
% of these variables is available in the documentation. 
% 
% The same inputs need to be set in the driver file StartMgPipe when running 
% mgPipe outside of this tutorial or directly in the "initMgPipe" function.
% 
% OPTIONAL INPUTS

% path to csv file for stratification criteria (if empty or not existent no criteria is used)
infoFilePath = '';

% path to a model of the host (e.g., human) to be joined with the
% microbioomes (default: no host)
hostPath = '';

% name of objective function of organisms, default='EX_biomass(e)'
objre = 'EX_biomass(e)';

% strategy used to build personalzied models. If true: create a global
% setup model that is pruned, if false, create each personalized mdoel one
% by one. The latter is recommended if there are 500 or more individual
% organisms.
buildSetupAll = true;

% if to save models with diet constrains implemented (default=false)
saveConstrModels = false;

% number of cores dedicated for parallelization (default=2)
numWorkers = 4;

% to enable also rich diet simulations (default=false)
rDiet = false;

% to enable personalized diet simulations (default=false)
pDiet = false;

% the type of FVA function to use to solve (true=fastFVA,
% false=fluxVariability)
fvaType = true;

% to manually set the lower bound on flux through the community biomass
% reaction (default=0.4 mmol/person/day)
lowerBMBound = 0.4;

% to set whether existing simulation results are rewritten (default=false)
repeatSim = false;

% to set if the input medium should be adapted through the adaptVMHDietToAGORA
% function or used as is (default=true)                  
adaptMedium = true; 

%% Pipeline run
% Calling the function initMgPipe will execute Part 1 to 3 of the pipeline.

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, 'resPath', resPath, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'hostPath', hostPath, 'objre', objre, 'buildSetupAll', buildSetupAll, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers, 'rDiet', rDiet, 'pDiet', pDiet, 'fvaType', fvaType, 'lowerBMBound', lowerBMBound, 'repeatSim', repeatSim, 'adaptMedium', adaptMedium);

%% Computed outputs
% # *Metabolic diversity* The number of mapped organisms for each individual 
% compared to the total number of unique reactions (extrapolated by the number 
% of reactions of each organism).Please, note that bigger circles with a number 
% inside represent overlapping individuals for metabolic diversity. 
% # *Classical multidimensional scaling of each individual reactions repertoire*
% 
% Flux Variability analysis for all the exchange reactions of the diet and fecal 
% compartment is also computed and saved in a file called "simRes". Specifically 
% what computed and saved are:
% 
% # *ID* a vector containing the names of metabolites for which FVA of exchange 
% reactions was computed
% # *fvaCt* a cell array containing min flux trough uptake and max trough secretion 
% exchanges (later used for computing NMPCs)
% # *nsCT* a cell array containing max flux trough uptake and min trough secretion 
% exchanges
% # *presol* an array containing the value of objectives for each microbiota 
% model with rich and selected diet
% # *inFesMat* cell array containing the names of the microbiota models that 
% reported an infeasible status when solved for their objective 
% 
% Finally, the net uptake and secretion potential are computed in a metabolite 
% resolved manner and saved in the files 'netSecretionFluxes.csv' and 'netUptakeFluxes.csv'
% results folder. They  indicate the maximal uptakr and production, respectively, of each 
% metabolite and are computed as the absolute value of the sum of the maximal secretion 
% flux with the maximal uptake flux. The similarity of metabolic secretion profiles 
% (using the net secretion potential as features) between individuals is also 
% evaluated with classical multidimensional scaling. 
% 
%% Stratification of samples
% If metadata for the analyzed samples is available (e.g.,
% disease state), the samples can be stratified based on this
% classification. To provide metadata, prepare an input file as in the
% example 'sampInfo.csv'. The path to the file with sample information
% needs to be provided as the variable infoFilePath. Note that the group 
% classification into groups in sampInfo.csv is arbitrary and not
% biological meaningful.

infoFilePath='sampInfo.csv'; 

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, 'resPath', resPath, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'hostPath', hostPath, 'objre', objre, 'buildSetupAll', buildSetupAll, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers, 'rDiet', rDiet, 'pDiet', pDiet, 'fvaType', fvaType, 'lowerBMBound', lowerBMBound, 'repeatSim', repeatSim, 'adaptMedium', adaptMedium);

%% Statistical analysis and plotting of generated fluxes
% If sample information as in sampInfo.csv is provided (e.g., healthy vs.
% disease state), statistical analysis can be performed to identify whether
% net secretion fluxes, net uptake fluxes, and reaction abundances can be
% performed. If the analyzed samples can be divided into two groups,
% Wilcoxon rank sum test will be used. If there are three or more groups,
% Kruskal-Wallis test will be performed.
% Note that no statistically significant differences will be found in the
% tutorial example due to the small sample sizes.
% Moreover, violin plots that show the computed fluxes separate by group
% will be generated.

% Define the input variables.
% Path to file with sample information (required)
infoFilePath='sampInfo.csv';

% Header in the file with sample information with the stratification to 
% analyze (if not provided, the second column will be chosen by default)
sampleGroupHeaders={'Group'};
% sampleGroupHeaders can contain more than one entry if multiple columns 
% with sample information (e.g., disease state, age group) should be analyzed.

% path with results of mgPipe that will be analyzed
resPath = [tutorialPath filesep 'Results'];

% define where results will be saved (optional, default folders will be
% generated otherwise)
statPath = [tutorialPath filesep 'Statistics'];
violinPath = [tutorialPath filesep 'ViolinPlots'];

%% perform the statistical analysis and save the results
analyzeMgPipeResults(infoFilePath,resPath,'statPath', statPath, 'violinPath', violinPath, 'sampleGroupHeaders', sampleGroupHeaders);

% Afterwards, the results of the statistical analysis will be available in
% the folder "Statistics". The files ending in "Statistics.txt" contain the
% calculated p-values and test results. If there are any fluxes or
% abundances that significantly differed between groups, there will be
% files ending in "significantFeatures.txt" listing only these instances.
% Created violin plots will be found in the folder "ViolinPlots". There
% will be an image in png and pdf format for each predicted metabolite's
% uptake and secretion potential. There will also be a file ending in
% "All_plots.pdf" containing all plots.

%% Analysis of strain-level contributions
% For metabolites of particular interest (e.g., for which the
% community-wide secretion potential was significantly different between
% disease cases and controls), the strains consuming and secreting the
% metabolites in each sample may be computed. This will yield insight into
% the contribution of each strain to each metabolite. Note that for
% metabolites for which the community wide secretion potential did not
% differ, the strains contributing to metabolites may still be
% significantly different.

% We will define a list of metabolites to analyze. As an example, we will
% take acetate and formate.
metList = {'ac','for'};

% path with microbiome models generated through mgPipe that will be analyzed
mmPath = [tutorialPath filesep 'Results'];

% create a new folder where strain contributions will be saved
mkdir([tutorialPath filesep 'StrainContributions']);
contPath = [tutorialPath filesep 'StrainContributions'];

[minFluxes,maxFluxes,fluxSpans] = predictMicrobeContributions(mmPath, 'resPath', contPath, 'dietFilePath', dietFilePath, 'metList', metList, 'numWorkers', numWorkers, 'lowerBMBound', lowerBMBound, 'adaptMedium', adaptMedium);

% The output 'minFluxes' shows the fluxes in the reverse direction through
% all internal exchange reactions that had nonzero flux for each analyzed
% metabolite. The output 'maxFluxes' shows the corresponding forward
% fluxes. 'fluxSpans' shows the distance between minimal and maximal
% fluxes for each internal exchange reaction with nonzero flux for each
% metabolite.

% Afterwards, statistical analysis of the strain contributions can also be
% performed.
analyzeMgPipeResults(infoFilePath,contPath, 'statPath', statPath, 'violinPath', violinPath, 'sampleGroupHeaders', sampleGroupHeaders);

%% Computation of shadow prices
% Shadow prices are routinely retrieved with each flux balance analysis
% solution. Briefly,the shadow price is a measurement for the value of a
% metabolite towards the optimized objective function, which indicates 
% whether the flux through the objective function would increase or 
% decrease when the availability of this metabolite would increase by one 
% unit (Palsson B. Systems biology : properties of reconstructed networks).
% For microbiome community models created through mgPipe, this will enable
% us to determine which metabolites are bottlenecks for the community's
% potential to secrete a metabolite of interest. This was performed for
% bile acids in Heinken et al., Microbiome (2019) 7:75.

% We will compute the shadow prices for acetate and formate as an example.
objectiveList={
    'EX_ac[fe]'
    'EX_for[fe]'
    };

% Here, we will compute all shadow prices that are nonzero. Thus, this
% include both metabolites that would increase and decrease the flux
% through the objective function if their availability was increased.
% Note that the definition of the shadow price depends on the solver. 
% To check the shadow price definitions for each solver, run the test
% script testDualRCostDefinition.
SPDef = 'Nonzero';

% create a new folder where shadow prices will be saved
mkdir([tutorialPath filesep 'ShadowPrices']);
spPath = [tutorialPath filesep 'ShadowPrices'];

[objectives,shadowPrices]=analyseObjectiveShadowPrices(mmPath, objectiveList, 'resultsFolder', spPath, 'SPDef', SPDef, 'numWorkers', numWorkers, 'dietFilePath', dietFilePath, 'lowerBMBound', lowerBMBound, 'adaptMedium', adaptMedium);

% Similar to the previous results, we can also perform statistical analysis
% on the computed shadow prices.
analyzeMgPipeResults(infoFilePath,spPath,'statPath', statPath, 'violinPath', violinPath,'sampleGroupHeaders', sampleGroupHeaders);

%% Using mgPipe to model host-microbiome co-metabolism
% If desired, the sample-specific microbiome models created by mgPipe can
% also be joined with a model of the host, e.g., with the human
% reconstruction Recon3D. 
% Let us create the same four microbiome models as before joined with
% Recon3D.
%
system('curl -LJO https://www.vmh.life/files/reconstructions/Recon/3D.01/Recon3D_301.zip')
unzip('Recon3D_301')
hostPath = [pwd filesep 'Recon3D_301' filesep 'Recon3DModel_301.mat'];
%
% Since host metabolites can now enter from the host model itself, the 
% adaptMedium input can be set to false.                 
adaptMedium = false; 
%
% overwrite the previous simulation for the purpose of the tutorial
repeatSim = true;
%
% If a host model is entered, it is also highly recommended to enter the
% host biomass reaction to generate coupling constraints for the host.
hostBiomassRxn = 'biomass_reaction';

% Run the pipeline including the host. Note that this will be more
% computationally intensive and take some time.

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, 'resPath', resPath, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'hostPath', hostPath, 'hostBiomassRxn', hostBiomassRxn, 'objre', objre, 'figForm', figForm, 'numWorkers', numWorkers, 'autoFix', autoFix, 'compMod', compMod, 'rDiet', rDiet, 'pDiet', pDiet, 'extSolve', extSolve, 'fvaType', fvaType, 'lowerBMBound', lowerBMBound, 'repeatSim', repeatSim, 'adaptMedium', adaptMedium);

