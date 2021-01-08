%% Refinement of genome-scale reconstructions through the DEMETER pipeline
%% Author: Almut Heinken, PhD, National University of Ireland Galway
% 01/2021

%% Introduction
% This tutorial introduces the semi-automated reconstruction pipeline  
% DEMETER (Data-drivEn METabolic nEtwork Refinement)  and its use for
% semi-automated, data-driven reconstruction of microbial strains.
% DEMETER was first applied to the reconstruction of 773 human gut 
% microbial organisms [1]. More recently, DEMETER was used to reconstruct
% 7,206 human microbial strains, resulting in AGORA2, an update to AGORA
% both in size and scope [2].
% The tutorial demonstrates how to use DEMETER to build a refined 
% genome-scale reconstruction for one or more target organisms from a draft 
% reconstruction retrieved from KBase.
% Steps include input data collection, input data integration, running the
% refinement pipeline, and testing the refined reconstruction against the
% input data. Steps Figure 1 summarizes the steps used to create a curated
% reconstruction.

%%

% Figure 1: Overview of the DEMETER pipeline [2]. The pipeline consists of 
% # collection of draft reconstructions, comparative genomic data, 
% biochemical and physiological data, and drug structures and microbial 
% conversion reactions, # conversion of data into a MATLAB-readable 
% format, and integration into pipeline functions, # data-driven 
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
% the apps at kbase.us/. A template narrative to facilitate this is 
% available at xx. 
% Next, please place all draft reconstructions that you want to refine into 
% one folder.
% For the sake of this tutorial, we will use two AGORA2 strains as well as
% six example strains that were randomly chosen from genomes present at 
% KBase and are not present in AGORA2. Draft reconstructions for these 
% strains have already been retrieved from KBase and placed into the folder 
% cobratoolbox/papers/2021_demeter/exampleDraftReconstructions.
% Alternatively, they are available in the template narrative xx.

% To start DEMETER, define the path to the folder where the draft 
% reconstructions are located.

initCobraToolbox
global CBTDIR
draftFolder = [CBTDIR filesep 'papers' filesep '2021_demeter' filesep 'exampleDraftReconstructions'];

%% Taxonomic and strain information
% Information on the taxonomic classification and organism properties
% (e.g., gram status, oxygen status, metabolism type, pathogenicity) has
% been collected for all 7,206 AGORA2 strains and is available in the file
% "input/AGORA2_infoFile.xlsx". It is highly recommmended to collect
% taxonomic information for your organism as this enables the propagation
% of experimental data from related strains that are already present in
% AGORA2. Taxonomic information can be retrieved from NCBI Taxonomy:
% https://www.ncbi.nlm.nih.gov/taxonomy/.
% Note that this step can be skipped, but will prevent informing the
% refined reconstructions through available experimental data.
% It is also recommended to find out whether your organism is gram-positive
% or -negative as this will enable curation of the biomass objective
% function. Gram status can be retrieved from literature searches or by
% querying the Integrated Microbial Genomes Database: 
% https://img.jgi.doe.gov/.
% To provide taxonomic information for your organisms, prepare a file in
% the same format as
% cobratoolbox/papers/2021_demeter/example_infoFile.xlsx.
% Please provide the path to the file with taxonomic information as the
% variable infoFilePath.
infoFilePath = [CBTDIR filesep 'papers' filesep '2021_demeter' filesep 'example_infoFile.xlsx'];

%% Step 2: Data integration
%% Propagating existing experimental data
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
% data afterwards as needed. The input file with taxonomic information will
% also be supplemented with gram staining information if possible, which
% serves to additionally inform and refine the reconstructions.
%% Using compatative genomic data as input data for DEMETER
% It is possible to inform the reconstruction with comparative
% genomics data retrieved from PubSEED. This will require either retrieving
% publicly available spreadsheets from PubSEED or the generation of such
% spreadsheets through PubSEED (https://pubseed.theseed.org). See 
% cobratoolbox/papers/2021_demeter/exampleSpreadsheets for examples of the
% required format.

% Annotated genes need to be provided as columns and mapped to
% corresponding reactions through the file
% cobratoolbox/papers/2021_demeter/input/InReactions.txt
% If comparative genomics data is available for your organism(s), please
% place it in a folder and provide the path to the folder as the variable
% spreadsheetFolder.
spreadsheetFolder = [CBTDIR filesep 'papers' filesep '2021_demeter' filesep 'exampleSpreadsheets'];

% Additionally, the PubSEED IDs for each strain to refine, as well as the 
% annotation version, need to be added to the "PubSeedID" and
% "Annotation version ID" columns, respectively, in the taxonomic 
% information file as shown in example_infoFile.xlsx.

% Define the path to a folder where the propagated data that will serve as 
% input for DEMETER should be stored (optional). Otherwise, a default path
% will be used.
inputDataFolder = [pwd filesep 'InputData'];

% To propagate available gram staining information and experimental data to 
% your organism(s), use the code
[adaptedInfoFilePath,inputDataFolder] = prepareInputData(infoFilePath,'inputDataFolder', inputDataFolder, 'spreadsheetFolder',spreadsheetFolder);

% The variable adaptedInfoFilePath contains the path to the taxonomic
% information file adapted with gram staining information.

%% Manual gathering of experimental data as input data for DEMETER
% In the folder inputDataFolder, you will find experimental data for your
% organism(s) retrieved and propagated from the already available data that
% had been collected for the AGORA2 resource. It is highly recommended that
% you manually inspect the propagated experimental data for your
% organism(s) at this point to verify it as organisms may differ on the 
% strain level. Moreover, it is recommended that you perform additional
% literature searches to inform the refinement of your organism(s).
% As you can see from inspection of the newly generated text files,
% for three of the four new strains used in this tutorial, experimental 
% data could be propagated from already reconstructed strains of the same
% species. For the fourth strain, Acetitomaculum ruminis DSM 5522, no
% related strain in present in AGORA2. To find experimental data on this
% organism, we go to PubMed (https://pubmed.ncbi.nlm.nih.gov/) and enter
% "Acetitomaculum ruminis" as a search term. This yields a published
% reference describing the metabolism of this species: PubMed ID 2500921,
% https://link.springer.com/article/10.1007/BF00416597.
% The paper reports that A. ruminis generates acetate through the acetogen
% pathway and uses formate and glucose as carbon sources. By converting
% this data into a MATLAB-readable format, we can ensure the associated 
% pathways are present in the refined reconstruction.

% To incorporate this data into the files with experimental data, use the 
% following code.

% read the file with carbon source information
data=readtable([inputDataFolder filesep 'CarbonSourcesTable.txt'], 'Delimiter', 'tab', 'ReadVariableNames', false);
data=table2cell(data);

% add information that glucose is used as a carbon source
findRow=find(strcmp(data(:,1),'Acetitomaculum_ruminis_DSM_5522'));
findCol=find(strcmp(data(1,:),'D-glucose'));
data{findRow,findCol}='1';
writetable(cell2table(data),[inputDataFolder filesep 'CarbonSourcesTable'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');

% read the file with information on consumed metabolites
data=readtable([inputDataFolder filesep 'uptakeTable.txt'], 'Delimiter', 'tab', 'ReadVariableNames', false);
data=table2cell(data);

% add information that formate is consumed
findRow=find(strcmp(data(:,1),'Acetitomaculum_ruminis_DSM_5522'));
findCol=strcmp(data(1,:),'Formate');
data{findRow,findCol}='1';
writetable(cell2table(data),[inputDataFolder filesep 'uptakeTable'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');

% read the file with fermentation product information
data=readtable([inputDataFolder filesep 'FermentationTable.txt'], 'Delimiter', 'tab', 'ReadVariableNames', false);
data=table2cell(data);

% add information that the acetogen pathway is used
findRow=find(strcmp(data(:,1),'Acetitomaculum_ruminis_DSM_5522'));
findCol=find(strcmp(data(1,:),'Acetogen pathway'));
data{findRow,findCol}='1';
writetable(cell2table(data),[inputDataFolder filesep 'FermentationTable'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');

% Alternatively, you can inspect the files manually and edit as needed.

%% Step 3: Iterative refinement 
%% Running the DEMETER pipeline
% With the experimental and comparative genomic data in place, the pipeline
% can be run. This will result in a refined reconstruction for each strain
% in the draft reconstruction folder as a mat file and optionally as a SBML
% file. Moreover, a summary of refinement steps, such as addition and
% removal of reactions during gap-filling, will be provided.
% Steps that are automatically carried out during this step include #
% translation of KBase reaction and metabolite to Virtual Metabolic Human 
% (https://www.vmh.life) nomenclature, # expansion of the reconstruction
% by pathways supported by experimental data, # refinement of pathways and
% gene rules against comparative genomic data, and # quality assurance
% and quality control ensuring e.g., thermodynamic feasibility.

% First, we will define the folders where pipeline results will be stored
% (optional, otherwise, default paths will be used).

% Define the path to a folder where the refined reconstructions will be
% stored in mat file format
refinedFolder = [pwd filesep 'RefinedReconstructions'];

% Define the path to a folder where the refined reconstructions will be
% stored in SBML file format
sbmlFolder = [pwd filesep 'refinedReconstructions_SBML'];

% Define the path where translated versions of the draft reconstructions
% will be stored. These reconstructions will not undergo the pipeline
% except for the translation step, which will enable them to pass the
% DEMETER test suite.
translatedDraftsFolder = [pwd filesep 'translatedDraftsFolder'];

% Define the folder where reports on performed gap-filling, QA/QC,
% and expansion the reconstructions will be saved.
summaryFolder = [pwd filesep 'refinementSummary'];

% Define whether the refined reconstructions will additionally be exported
% as SBML files (default = false).
createSBML = true;

% Define the number of workers for parallel computing.
numWorkers = 4;

% Define a name for the reconstruction resource (optional)
reconVersion = 'TutorialExample';

% Run the pipeline.
runPipeline(draftFolder, 'infoFilePath', adaptedInfoFilePath, 'inputDataFolder', inputDataFolder, 'refinedFolder', refinedFolder, 'translatedDraftsFolder', translatedDraftsFolder, 'summaryFolder', summaryFolder, 'numWorkers', numWorkers, 'reconVersion', reconVersion, 'createSBML', createSBML, 'sbmlFolder', sbmlFolder)

% run test suite
% draft folder for testing
createMatfileKBaseDraftModels
draftFolder = [rootDir filesep 'Current_Version_AGORA2' filesep 'Draft_Reconstructions_Matfiles'];
runTestSuiteTools(draftFolder, refinedFolder, 'numWorkers', numWorkers, 'testResultsFolder', testResultsFolder, 'reconVersion', reconVersion)

% get model properties
draftFolder = [rootDir filesep 'Current_Version_AGORA2' filesep 'Draft_Reconstructions_Matfiles'];
computeModelProperties(draftFolder, refinedFolder, 'numWorkers', numWorkers, 'propertiesFolder', propertiesFolder, 'reconVersion', reconVersion)


