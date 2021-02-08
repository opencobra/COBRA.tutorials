%% Creation, interrogation, and analysis of personalized microbiota models from the AGORA models through metagenomic data integration
%% Author: Almut Heinken, National University of Ireland Galway and Federico Baldini, Luxembourg Centre for Systems Biomedicine 
%% INTRODUCTION
% This tutorial demonstrates the various function of the Microbiome Modeling 
% Toolbox for the creation, interrogation, analysis, and stratification of personalized 
% microbiome models. The tutorial consists of the following steps:
%% 
% # Mapping files with abundances obtained from metagenomic sequencing data 
% to the AGORA models
% # Creating pan-models from AGORA on different taxonomical levels
% # Building personalized models
% # Computing net fecal uptake and secretion potential for all microbial metabolites
% # Calculating subsystem abundances on the microbiome level
% # Clustering and statistical analysis of net uptake and net secretion fluxes 
% based on subgroups in samples
% # Targeted analysis: computation of strain-level contributions in each sample 
% to specific metabolites
% # Targeted analysis: computation of shadow prices that reveal metabolite dependencies 
% for a specific objective
% # Integration of a host organism with personalized microbiome models
%% 
% Please note that this tutorial uses as an example a small dataset (4 columns 
% and 10 rows) with the purpose of demonstrating the functionalities of the pipeline. 
% We recommend using high-performance computing clusters when assembling and simulating 
% from bigger datasets.
%% DRIVER
% The  pipeline driver script cobratoolbox/papers/2018_microbiomeModelingToolbox/startMgPipe.m 
% contains the neccessary inputs to run mgPipe for your dataset. Replace the input 
% variables in startMgPipe with your own dataset and modifiy as needed.
%% Initialize the COBRA Toolbox

initCobraToolbox(false)
%% 
% Set a solver (recommended: ibm_cplex)

solverOK=changeCobraSolver('ibm_cplex','LP');
%% Prepare input data and models
% We first set the paths to input and output files change directory to where 
% the tutorial is located

tutorialPath = fileparts(which('tutorial_mgPipe'));
cd(tutorialPath);
%% 
% We will use the AGORA resource (Magnusdottir et al., Nat Biotechnol. 2017 
% Jan;35(1):81-89) in this tutorial. AGORA version 1.03 is available at <https://github.com/VirtualMetabolicHuman/AGORA 
% https://github.com/VirtualMetabolicHuman/AGORA>. Download AGORA and place the 
% models into a folder.

system('curl -LJO https://github.com/VirtualMetabolicHuman/AGORA/archive/master.zip')
unzip('AGORA-master')
modPath = [pwd filesep 'AGORA-master' filesep 'CurrentVersion' filesep 'AGORA_1_03' filesep' 'AGORA_1_03_mat'];
%% Preparing the input file with normalized abundances
% Preparing an input file suitable for mgPipe requires mapping the organisms 
% to the nomenclature of the AGORA models. A function that can help with this 
% is available in the COBRA Toolbox.
% 
% For an example of how to map metagenomic sequences to AGORA models, enter 
% the following code:

system('curl -LJO https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00001248/pipelines/3.0/file/SRP065497_taxonomy_abundances_v3.0.tsv')

[translatedAbundances,normalizedAbundances,unmappedRows] = translateMetagenome2AGORA('SRP065497_taxonomy_abundances_v3.0.tsv','Species');
%% 
% _translateMetagenome2AGORA_ translates the output of common sequencing pipelines 
% (e.g., MetaPhlAn) to the names of the AGORA models. However, the function may 
% not work for the nomenclature of some input file types it has never been tested 
% on. Manual inspection to map taxa may still be needed. Inspect the output variable 
% unmappedRows to check for mismatches in nomenclature. Such mismatches can be 
% correcting by adding nomenclature corrections after line 115 of the _translateMetagenome2AGORA_ 
% function. It is requested to also submit such nomenclature corrections as a 
% pull request.
% 
% If there are entries in unmappedRows because the strain or species is not 
% captured by the AGORA resource, you are directed to the tutorial for DEMETER 
% (cobratoolbox/tutorials/reconstruction/demeter). DEMETER is a efficient reconstruction 
% refinement pipeline that was previously used to reconstruct AGORA.  DEMETER 
% enables fast, data-driven reconstruction of new strains following the same quality 
% standards as AGORA.
% 
% To normalize the relative abundances from the input file that had been mapped 
% to AGORA organisms, run the function normalizeCoverage. It is recommended to 
% set the input variable cutoff (e.g., 0.0001), which removes all strains below 
% a certain normalized relative abundance in a sample. Removing taxa with very 
% low normalized coverage will result in faster simulations.

cutoff = 0.0001;
[normalizedCoverage,abunFilePath] = normalizeCoverage(translatedAbundances,cutoff);
%% 
% The file normalizedCoverage.csv now contains relative normalized abundances 
% mapped to the nomenclature of AGORA that could serve as input for mgPipe. The 
% output variable abunFilePath can directly serve as input for the corresponding 
% input variable in initMgPipe. 
%% Creating pan-models from AGORA on different taxonomical levels
% As the abundances in this example are provided on the species level, the AGORA 
% reconstructions, which are strain-specific, cannot serve as direct input. To 
% enable the creation of personalized microbiome models from species-level relative 
% abundances (or alternatively genus- or family-level), pan-models need to be 
% created from AGORA. The following code performs this on the species level. To 
% create pan-models on other species level up to phylum, modify the input variable 
% taxonLevel.

panPath=[tutorialPath filesep 'panSpeciesModels'];
mkdir(panPath)

taxonLevel='Species';

createPanModels(modPath,panPath,taxonLevel);
%% 
% By setting panPath as the input variable modPath for initMgPipe, personalized 
% microbiome models for the samples from the above study with the EMBL-EBI accession 
% ID SRP065497 could now readily  be created and interrogated through mgPipe. 
% As this would be too computationally intensive for the purpose of this tutorial, 
% this will not be performed here.
%% Building personalized models 
% While personalized models are built by executing the function _initMgPipe_, 
% the following steps are performed automatically.
%% 
% * Classicial multidimensional scaling (PCoA) of the relative abundances is 
% performed.
% * The diversity in terms of organisms present ine ach sample is plotted.
% * The relatives abundances of each reaction present in at least one organism 
% in each sample is computed (value of 0 to 1 with 1 indicating that every organism 
% in the sample has the reaction).
%% 
% First, we will set the path where to save results.

mkdir('Results');
resPath = [tutorialPath filesep 'Results'];
%% 
% Then, we will define the simulated dietary regime that will be implemented 
% on the personalized models. Here, we will use an "Average European" diet that 
% is located in the folder
% 
% cobratoolbox/papers/2018_microbiomeModelingToolbox/input
% 
% The folder also contains various other simulated diets, e.g., Vegan, gluten-free.

dietFilePath='AverageEuropeanDiet';
%% 
% Next, we set the path and the name of the file from which to load the abundances. 
% For this tutorial, to reduce the time of computations, we will use a reduced 
% version of the example file (normCoverageReduced.csv) provided in the folder 
% Resources: only 4 individuals and 10 strains will be considered.

abunFilePath='normCoverageReduced.csv';
%% 
% For the creation of personalized models, mgPipe offers two strategies that 
% are selected in the input variable buildSetupAll. By default, a global metabolic 
% model (setup) is constructed containing all the microbes listed in the study. 
% Personalized models are built for each sample by pruning the global model. As 
% another option, which is recommended for studies with 300 or more microbes total, 
% personalized models can be built serially. In both cases, a community biomass 
% reaction is built from the relative abundances in each sample.
% 
% By default, mgPipe will subsequently compute the metabolic profiles of each 
% personalized model. Briefly, the chosen dietary regime is implemented and the 
% total net uptake and net secretion potential for each metabolite that can be 
% transported by at least one microbe in the study is computed for each personalized 
% model. The PCoA of computed fluxes is also plotted. If this is not neccessary 
% and targeted analyses are instead planned (see later sections in the tutorial), 
% the computation can be disabled through the input variable computeProfiles. 
% In either case, the personalized models with dietary constraints implemented 
% can be exported through the variable saveConstrModels. Exporting the constrained 
% models will also enable customized analyses outside mgPipe.
% 
% To define whether flux variability analysis to compute the metabolic profiles 
% should be performed, and set the input variable computeProfiles. If fastFVA 
% is installed, it will be used, otherwise, fluxVariability will be used. It is 
% highly recommended to install fastFVA as this drastically reduces computation 
% times in large datasets.

computeProfiles = true;
%% Setting optional inputs
% Next inputs will define:
%% 
% # name of the objective function of organisms
% # whether to create the personalized models by pruning a global setup model 
% (default) or one by one
% # number of cores to use for the pipeline execution
% # if stratification criteria (e..g., cases and controls) are available for 
% the samples
% # if to simulate also a rich diet
% # if to save models with diet constraints implemented
% # the lower bound in mmol/person/day set on the community biomass reaction
% # whether already existing simulation results should be overwritten
% # whether the simulated dietary regime should be supplemented with compounds 
% typically found in the gut thus mimicking the human intestinal environment
%% 
% The following setting should work for almost any system, but please check 
% carefully to be sure these options are valid for you. A more detailed description 
% of these variables is available in the documentation. 
% 
% The same inputs need to be set in the driver file StartMgPipe when running 
% mgPipe outside of this tutorial or directly in the "initMgPipe" function.
% 
% OPTIONAL INPUTS
% 
% path to csv file for stratification criteria (if empty or not existent no 
% criteria is used)

infoFilePath = '';
%% 
% path to a model of the host (e.g., human) to be joined with the microbiomes 
% (default: no host)

hostPath = '';
%% 
% name of objective function of organisms, default='EX_biomass(e)'

objre = 'EX_biomass(e)';
%% 
% strategy used to build personalized models. If true: create a global setup 
% model that is pruned, if false, create each personalized model one by one. The 
% latter is recommended if there are 300 or more individual organisms.

buildSetupAll = true;
%% 
% if to save models with diet constrains implemented (default=false)

saveConstrModels = true;
%% 
% number of cores dedicated for parallelization (default=2)

numWorkers = 4;
%% 
% to enable also rich diet simulations (default=false)

rDiet = false;
%% 
% to enable personalized diet simulations (default=false)

pDiet = false;
%% 
% to manually set the lower bound on flux through the community biomass reaction 
% (default=0.4 mmol/person/day)

lowerBMBound = 0.4;
%% 
% to set whether existing simulation results are rewritten (default=false)

repeatSim = false;
%% 
% to set if the input medium should be adapted through the adaptVMHDietToAGORA 
% function or used as is (default=true)                  

adaptMedium = true;
%% Pipeline run
% Calling the function initMgPipe will execute Part 1 to 3 of the pipeline.

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, computeProfiles, 'resPath', resPath, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'hostPath', hostPath, 'objre', objre, 'buildSetupAll', buildSetupAll, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers, 'rDiet', rDiet, 'pDiet', pDiet, 'lowerBMBound', lowerBMBound, 'repeatSim', repeatSim, 'adaptMedium', adaptMedium);
%% Computed outputs
%% 
% # *Metabolic diversity* The number of mapped organisms for each individual 
% compared to the total number of unique reactions (extrapolated by the number 
% of reactions of each organism).Please, note that bigger circles with a number 
% inside represent overlapping individuals for metabolic diversity. 
% # *Classical multidimensional scaling of each individual reactions repertoire*
%% 
% *If flux variability analysis is performed and net uptake and secretion potential 
% are computed:*
% 
% Flux Variability analysis for all the exchange reactions of the diet and fecal 
% compartment are performed and temporarily saved in a file called "simRes". Specifically 
% what is temporarily saved is:
%% 
% # *fvaCt* a cell array containing min flux through uptake and max trough secretion 
% exchanges
% # *nsCT* a cell array containing max flux through uptake and min trough secretion 
% exchanges
% # *presol* an array containing the value of objectives for each microbiota 
% model with rich and selected diet
% # *inFesMat* cell array containing the names of the microbiota models that 
% reported an infeasible status when solved for their objective 
%% 
% Finally, the net uptake and secretion potential are computed in a metabolite 
% resolved manner and saved in the output variables 'netSecretionFluxes' and netUptakeFluxes, 
% and the files 'netSecretionFluxes.csv' and 'netUptakeFluxes.csv' in the results 
% folder. They indicate the maximal uptake and production, respectively, of each 
% metabolite and are computed as the absolute value of the sum of the maximal 
% secretion flux with the maximal uptake flux. The similarity of metabolic secretion 
% profiles (using the net secretion potential as features) between individuals 
% is also evaluated with classical multidimensional scaling. 
% 
% In the export of models with dietary constraints is desired, they will be 
% found in the folder constrModelsPath.
%% Calculating subsystem abundances
% The out file "reactions.csv" in the Results folder cotnain the relative abundance 
% of each reaction in each sample. A description of each reaction can be found 
% by searching for the reaction ID in the file
% 
% cobratoolbox/papers/2021_demeter/input/ReactionDatabase.txt
% 
% or on the Virtual Metabolic Human website (https://www.vmh.life/).
% 
% For a more convenient overview of pathways altered in their relative abundances 
% across samples, the relative abundances can also be calculated for each sample 
% on the subsystem level by entering the code

reactionAbundancePath=[resPath filesep 'reactions.csv'];
subsystemAbundance = calculateSubsystemAbundance(reactionAbundancePath);
%% 
% The calculated subsystem abundances can be found in the output variable subsystemAbundance 
% and in the file "SubsystemAbundance.txt". A subsystem abundance of 1 indicates 
% that each reaction in the subsystem found in at least one sample is present 
% at a relative abundance of 1 in the sample.
%% Stratification of samples
% If metadata for the analyzed samples is available (e.g., disease state), the 
% samples can be stratified based on this classification. To provide metadata, 
% prepare an input file as in the example 'sampInfo.csv'. The path to the file 
% with sample information needs to be provided as the variable infoFilePath. Note 
% that the group classification into groups in sampInfo.csv is arbitrary and not 
% biological meaningful.
% 
% Run the pipeline agian, this time with sample stratification into groups provided. 
% The PCoAs will now be labelled with sample stratification into groups.

infoFilePath='sampInfo.csv'; 
saveConstrModels = false;

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, computeProfiles, 'resPath', resPath, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'hostPath', hostPath, 'objre', objre, 'buildSetupAll', buildSetupAll, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers, 'rDiet', rDiet, 'pDiet', pDiet, 'lowerBMBound', lowerBMBound, 'repeatSim', repeatSim, 'adaptMedium', adaptMedium);
%% Statistical analysis and plotting of generated fluxes
% If sample information as in sampInfo.csv is provided (e.g., healthy vs. disease 
% state), statistical analysis can be performed to identify whether net secretion 
% fluxes, net uptake fluxes, and reaction abundances can be performed. If the 
% analyzed samples can be divided into two groups, Wilcoxon rank sum test will 
% be used. If there are three or more groups, Kruskal-Wallis test will be performed. 
% Note that no statistically significant differences will be found in the tutorial 
% example due to the small sample sizes. Moreover, violin plots that show the 
% computed fluxes separate by group will be generated.
% 
% Define the input variables. Path to file with sample information (required)

infoFilePath='sampInfo.csv';
%% 
% Header in the file with sample information with the stratification to  analyze 
% (if not provided, the second column will be chosen by default)

sampleGroupHeaders={'Group'};
%% 
% sampleGroupHeaders can contain more than one entry if multiple columns  with 
% sample information (e.g., disease state, age group) should be analyzed.
% 
% path with results of mgPipe that will be analyzed

resPath = [tutorialPath filesep 'Results'];
%% 
% define where results will be saved (optional, default folders will be generated 
% otherwise)

statPath = [tutorialPath filesep 'Statistics'];
violinPath = [tutorialPath filesep 'ViolinPlots'];
%% 
% *To perform the statistical analysis and save the results, enter the code*

analyzeMgPipeResults(infoFilePath,resPath,'statPath', statPath, 'violinPath', violinPath, 'sampleGroupHeaders', sampleGroupHeaders);
%% 
% Afterwards, the results of the statistical analysis will be available in the 
% folder "Statistics". The files ending in "Statistics.txt" contain the calculated 
% p-values and test results. If there are any fluxes or abundances that significantly 
% differed between groups, there will be files ending in "significantFeatures.txt" 
% listing only these instances. Created violin plots will be found in the folder 
% "ViolinPlots". There will be an image in png and pdf format for each predicted 
% metabolite's uptake and secretion potential. There will also be a file ending 
% in "All_plots.pdf" containing all plots.
%% Targeted analysis: Strain-level contributions to metabolites of interest
% For metabolites of particular interest (e.g., for which the community-wide 
% secretion potential was significantly different between disease cases and controls), 
% the strains consuming and secreting the metabolites in each sample may be computed. 
% This will yield insight into the contribution of each strain to each metabolite. 
% Note that for metabolites for which the community wide secretion potential did 
% not differ, the strains contributing to metabolites may still be significantly 
% different.
% 
% The first step for the preparation of targeted analyses is the export of models 
% that had already been constrained with the simulated dietary regime. This had 
% already been done above through the saveConstrModels input in mgPipe above. 
% Now, will set the input variable modPath to the folder with personalized models 
% constrained with the simulated diet.

constrModPath = [resPath filesep 'Diet'];
%% 
% We will define a list of metabolites to analyze. As an example, we will take 
% acetate and formate.

metList = {'ac','for'};
%% 
% create a new folder where strain contributions will be saved

mkdir([tutorialPath filesep 'StrainContributions']);
contPath = [tutorialPath filesep 'StrainContributions'];

[minFluxes,maxFluxes,fluxSpans] = predictMicrobeContributions(constrModPath, 'resPath', contPath, 'metList', metList, 'numWorkers', numWorkers);
%% 
% The output 'minFluxes' shows the fluxes in the reverse direction through all 
% internal exchange reactions that had nonzero flux for each analyzed metabolite. 
% The output 'maxFluxes' shows the corresponding forward fluxes. 'fluxSpans' shows 
% the distance between minimal and maximal fluxes for each internal exchange reaction 
% with nonzero flux for each metabolite.
% 
% Afterwards, statistical analysis of the strain contributions can also be performed.

analyzeMgPipeResults(infoFilePath,contPath, 'statPath', statPath, 'violinPath', violinPath, 'sampleGroupHeaders', sampleGroupHeaders);
%% Targeted analysis: Computation of shadow prices for secreted metabolites of interest
% Shadow prices are routinely retrieved with each flux balance analysis solution. 
% Briefly, the shadow price is a measurement for the value of a metabolite towards 
% the optimized objective function, which indicates whether the flux through the 
% objective function would increase or decrease when the availability of this 
% metabolite would increase by one unit (Palsson B. Systems biology : properties 
% of reconstructed networks). For microbiome community models created through 
% mgPipe, this will enable us to determine which metabolites are bottlenecks for 
% the community's potential to secrete a metabolite of interest. This was performed 
% for bile acids in Heinken et al., Microbiome (2019) 7:75. However, any other 
% reaction in the personalized model (or any other model) can also serve as the 
% input.
% 
% We will compute the shadow prices for acetate and formate as an example.

objectiveList={
    'EX_ac[fe]'
    'EX_for[fe]'
    };
%% 
% Here, we will compute all shadow prices that are nonzero. Thus, this include 
% both metabolites that would increase and decrease the flux through the objective 
% function if their availability was increased. Note that the definition of the 
% shadow price depends on the solver.  To check the shadow price definitions for 
% each solver, run the test script testDualRCostDefinition.

SPDef = 'Nonzero';
%% 
% create a new folder where shadow prices will be saved

mkdir([tutorialPath filesep 'ShadowPrices']);
spPath = [tutorialPath filesep 'ShadowPrices'];

[objectives,shadowPrices]=analyseObjectiveShadowPrices(constrModPath, objectiveList, 'resultsFolder', spPath, 'SPDef', SPDef, 'numWorkers', numWorkers);
%% 
% Similar to the previous results, we can also perform statistical analysis 
% on the computed shadow prices.

analyzeMgPipeResults(infoFilePath,spPath,'statPath', statPath, 'violinPath', violinPath,'sampleGroupHeaders', sampleGroupHeaders);
%% Using mgPipe to model host-microbiome co-metabolism
% If desired, the sample-specific microbiome models created by mgPipe can also 
% be joined with a model of the host, e.g., with the human reconstruction Recon3D. 
% Let us create the same four microbiome models as before joined with Recon3D.

system('curl -LJO https://www.vmh.life/files/reconstructions/Recon/3D.01/Recon3D_301.zip')
unzip('Recon3D_301')
hostPath = [pwd filesep 'Recon3D_301' filesep 'Recon3DModel_301.mat'];
%% 
% First, we will delete the previously created setup model so it can be overwritten/

delete([resPath filesep 'mapInfo.mat'],[resPath filesep 'Setup_allbacs.mat'])
%% 
% Since host metabolites can now enter from the host model itself, the  adaptMedium 
% input can be set to false.                 

adaptMedium = false; 
%% 
% If a host model is entered, it is also highly recommended to enter the host 
% biomass reaction to generate coupling constraints for the host.

hostBiomassRxn = 'biomass_reaction';
%% 
% The upper bound on the flux through the host biomass reaction can also be 
% constrained by entering the input variable hostBiomassRxnFlux (default: 1).

hostBiomassRxnFlux = 1;
%% 
% For the sake of this tutorial, we will only create personalized models combined 
% with Recon3D and constrained with the Average European diet, not run FVA again. 
% Naturally, all analyses shown in this tutorial can also be performed with the 
% host present.

computeProfiles = false;
%% 
% Run the pipeline to create personalized host-microbiome models.

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, abunFilePath, computeProfiles, 'resPath', resPath, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'hostPath', hostPath, 'hostBiomassRxn', hostBiomassRxn, 'hostBiomassRxnFlux', hostBiomassRxnFlux, 'objre', objre, 'buildSetupAll', buildSetupAll, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers, 'rDiet', rDiet, 'pDiet', pDiet, 'lowerBMBound', lowerBMBound, 'repeatSim', repeatSim, 'adaptMedium', adaptMedium);