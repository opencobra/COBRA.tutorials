%% Creation of curated genome-scale reconstructions through a semi-automatic refinement pipeline
%% Almut Heinken, PhD
%% National University of Ireland Galway

%% Introduction
% This tutorial introduces the AGORA2 pipeline and its use for
% semi-automated, data-driven reconstruction of microbial strains as first
% introduced in [1]. The tutorial demonstrates how to build a refined 
% genome-scale reconstruction for one or more target organisms from a draft 
% reconstruction retrieved from KBase.
% Steps include input data collection, input data integration, running the
% refinement pipeline, and testing the refined reconstruction against the
% input data. Steps Figure 1 summarizes the steps used to create a curated
% reconstruction.

% Figure 1: Overview of the AGORA2 pipeline. The pipeline consists of (i) 
% collection of draft reconstructions, comparative genomic data, 
% biochemical and physiological data, and drug structures and microbial 
% conversion reactions, (ii) conversion of data into a MATLAB-readable 
% format, and integration into pipeline functions, (iii) data-driven 
% refinement of KBase draft reconstructions with continuous testing 
% against input data and subsequent update to pipeline functions to 
% correct for false predictions.

%% Step 1: Data collection
% This part explains how to collect input data and convert it into a format
% readible by the pipeline. The minimal input for building a refined
% reconstruction is a KBase draft reconstruction. Recommended inputs
% include taxonomic information on the target organism, experimental data
% if the species is not yet found in the already collected data, and
% comparative genomic analyses.
%% KBase draft reconstruction-required
% Create KBase draft reconstruction(s) for the target organism(s) by using
% the apps at kbase.us/. The process is further explained in the tutorial 
% narratives at KBase as well as in the manuscript version of this tutorial
% available at xx. Please place all draft reconstructions that you want to 
% refine into one folder.
% For the sake of this tutorial, we will use four example strains that were 
% randomly chosen from genomes present at KBase and are not present in 
% AGORA2. Draft reconstructions for these strains have already been
% retrieved from KBase and placed into the "examples" folder.
% Define the path to the folder where the draft reconstructions are
% located.
fileDir = fileparts(which('ReactionTranslationTable.txt'));
inputFolder = strrep(fileDir,'input','examples');

%% Taxonomic and strain information
% Information on the taxonomic classification and organism properties
% (e.g., gram status, oxygen status, metabolism type, pathogenicity) has
% been collected for all 7,206 AGORA2 strains and is available in the file
% "input/AGORA2_infoFile.xlsx". It is highly recommmended to collect
% taxonomic information for your organism as this enables the propagation
% of experimental data from related strains that are already present in
% AGORA2. Taxonomic information can be retrieved from NCBI Taxonomy:
% https://www.ncbi.nlm.nih.gov/taxonomy/.
% It is also recommended to find out whether your organism is gram-positive
% or -negative as this will enable curation of the biomass objective
% function. Gram status can be retrieved from literature searches or by
% querying the Integrated Microbial Genomes Database: 
% https://img.jgi.doe.gov/.
% To add information to "AGORA2_infoFile.xlsx", run the code

%% Propgagating existing experimental data
% The AGORA2 reconstruction pipeline contains experimental data for >1,500
% microbial species for (i) carbon sources, (ii) fermentation pathways,
% (iii) growth requirements, (iv) consumed metabolites, and (v) secreted
% metabolites. The data are stored in MATLAB-readable format in the input
% folder as the files CarbonSourcesTable.txt, FermentationTable.txt, 
% NutrientRequirementsTable.txt, secretionProductTable.txt, and
% uptakeTable.txt. If a strain related to your organism is already present
% in the AGORA2 resource and has data available, this data can be 
% propagated to your new organism. This requires taxonomic information 
% organized in a table as described in the previous step. Note that you can 
% still adapt the experimental data to strain-specific findings and add new 
% data afterwards as needed.
% To propagate available experimental data to your organism(s), use the
% code
[propagatedData] = propagateExperimentalData(inputData,strainInformation);

%% Gathering experimental data
% As you can see from inspection of the newly generated text files,
% for three of the four strains used in this tutorial, experimental data
% could be propagated from already reconstructed strains of the same
% species.For the fourth strain, Acetitomaculum ruminis DSM 5522, no
% related strain in present in AGORA2. To find experimental data on this
% organism, we go to PubMed (https://pubmed.ncbi.nlm.nih.gov/) and enter
% "Acetitomaculum ruminis" as a search term. This yields a published
% reference describing the metabolism of this species: PubMed ID 2500921,
% https://link.springer.com/article/10.1007/BF00416597.
% The paper tells us that A. ruminis generates acetate through the acetogen
% pathway and uses formate and glucose as carbon sources. By converting
% this data into a MATLAB-readable format, we can ensure the associated 
% pathways are present in the curated reconstruction.
% To incorporate this data into CarbonSourcesTable.txt and 
% FermentationTable.txt, respectively, use the code

%% Comparative genomic analyses

%
%
% Based on the nutrient requirement data for the microbe (table prepared
% above), determine whether the _in vitro_ data and _in silico_ simulations of
% the microbe match (nutrient essential or non-essential).
%
% IMPORTANT: It is necessary to check the model constraints carefully before
% performing this step. If there are any reactions that should or should not always
% carry a non-zero flux, then such constraints need to be applied or removed before
% running the nutrient essentiality test since such constraints can affect the
% _in silico_ nutrient essentiality results.
%
% Collect your data and prepare the data files as explained in the tutorial
% "PrepareInputData".
%
% _*IMPORTANT*: _Make sure to run the function "|validateExperimentalInputFiles|"
% on your prepared files to make sure they have the right format for the pipeline
% to run without errors.
% validateExperimentalInputFiles

% home = getenv('HOME');
% w = warning ('off','all');

%% initialize the COBRA Toolbox and solvers
initCobraToolbox
% solverOK=changeCobraSolver('gurobi','LP');
solverOK=changeCobraSolver('ibm_cplex','LP');
% prevent creation of log files
changeCobraSolverParams('LP', 'logFile', 0);



%% set the input and output folders
inputFolder=[rootDir filesep 'Current_Version_AGORA2' filesep 'Input_Models' filesep];

%% New run %%
outputFolder=[rootDir filesep 'Current_Version_AGORA2' filesep 'Output_Models' filesep];

numWorkers=8;

%% start the pipeline
pipelineDriver

