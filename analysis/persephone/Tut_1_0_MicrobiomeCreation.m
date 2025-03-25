% Tutorial 1: creating Human-microbiome models
% Authors: Bram Nap -  07-2024
% Welcome to the first tutorial in this three series tutorial set on how to create and analyse whole-body models (WBMs) and host-microbiome whole body models (mWBMs). In this tutorial we will cover how to create mWMBs from metagenomic reads data. We will explain how to use Microbial Abundances Retrieved from Sequencing data—automated NCBI Taxonomy (MARS) [1] to map metagenomic reads to AGORA2 database [2]. We will then explain how to process the output of MARS through the Microbiome Modelling Toolbox [3] to generate community microbiome models. Finally we will explain how to link the community microbiome models with the whole-body models [4] and ensure their feasibility.
% 
% Section 1: Setup
% Here we show which steps are required to set up your device to support the various functions used to create and analyse human-microbiome models for all three tutorials. First we need to have a copy of COBRA toolbox. COBRA Toolbox download instructions can be found at https://github.com/opencobra/cobratoolbox
% To see if the COBRA toolbox is installed correctly we run initCobraToolbox
global CBTDIR
if isempty(CBTDIR)
    initCobraToolbox
end
% As the function gives no errors or warnings we know that the COBRA toolbox is correctly set-up and ready to use.
% 
% To decrease simulation times of the models we need an industrial solver. The different solvers supported are
% ibm_cplex
% tomlab_cplex
% gurobi
% mosek
% To see MATLAB version and solver version compatibility see https://opencobra.github.io/cobratoolbox/stable/installation.html#solver-installation. Various solvers can be obtained through an free academic license. For modelling purposes here, we recommomend using ibm_cplex if possible.
% To set our solver we use changeCobraSolver. If another solver than ibm_cplex is used, replace the 'ibm_cplex' in the the code line below with the solver name you use as found in the bullet list above.
solver = 'ibm_cplex';

changeCobraSolver(solver);

% MATLAB also needs to have two additional toolboxes installed, the Statistics and Machine Learning Toolbox and the Parallel Computing Toolbox. We can check the installations with the code matlab.addons.isAddonEnabled.
% Check the install of the parallel computing toolbox
parallelToolboxInstall = matlab.addons.isAddonEnabled('Parallel Computing Toolbox')

% Check the install of the statistics and machine learning toolbox
statToolboxInstall = matlab.addons.isAddonEnabled('Statistics and Machine Learning Toolbox')
% If the parallelEnabled and statisticsEnables are both 1 (true) both toolboxes are correctly installed in MATLAB and ready to be used. If either one is 0 we recommended adding the toolbox from Home tab -> Add-Ons -> Get Add-Ons. The parallel computing toolbox is a must as the code creating the mWMBs does not work without out. The statistics toolbox is not required to create the mWBMs but without it we cannot run the statistical analyses performed in tutorial 2.
% Before we can start with creating mWBMs, we need to set up the paths to the general results directory, the metadata file, and the taxonomy assigned reads table. We will use the general results directory to create a folder structure where we can easily store our results in pre-defined locations. This will help as functions generally need the outputs of functions run before them. The metadata and reads tables undergo sanity checks. The variables we will set up are then:
% resultDir - The path where all the results of this tutorial will be stored. This ensure all your results will be in one place and are thus easily accesible. Make sure the the folder is accesible to you.
% metadataPath - The location of the metadata file. This is used in both the generation of human-microbiome models as well as the analysis of results towards the end of the tutorial.
% readsTablePath - The path to the the taxonomic reads file containing the amount of reads per taxonomic assignment. If the taxonomic assignments are in the reads file the first column needs to be called Taxon. If not the first column of the readsFile and the first column of taxaTable needs to have to same name. See the example files.
% We set up the variables in the paths variable, which is a structure. Each field in the structure is dedicated to store the relevant variables used in the different steps of Persephone. Normally we would define each variable at the start in the configuration file. However here we will only specifiy the output directories used in Persephone. The other variables we will define as we go along for clarity.
% 
% Please copy and paste the paths to the required files in the code below. Remember to add file extensions where relevant.
%Define the location of the required files, make sure you dont forget the file extensions where relevant!
resultPath = 'D:\OneDrive - National University of Ireland, Galway\Protocal_Pipeline\Test20252701';

paths.General.metadataPath = 'D:\OneDrive - National University of Ireland, Galway\Protocal_Pipeline\InputsTutorial\parkinsonSample5\metadataSubset5.xlsx';

paths.Mars.readsTablePath = 'D:\OneDrive - National University of Ireland, Galway\Protocal_Pipeline\InputsTutorial\parkinsonSample5\microbiomeSubset5.xlsx';
% Now we will define some paths as to where we want to store all of our results. We do not have to make these directories ourselves, the code will do it for us.
paths.seqC.outputPathSeqC = [resultPath, filesep, 'resultSeqC'];
paths.Mars.outputPathMARS = [resultPath, filesep, 'resultMars'];
paths.mgPipe.mgpipePath = [resultPath, filesep, 'resultMgPipe'];
paths.persWBM.personalizeWBMPath = [resultPath, filesep, 'personalisedWBMs'];
paths.mWBM.mWBMPath = [resultPath, filesep, 'mWBMmodels'];
paths.fba.fluxPath = [resultPath, filesep, 'resultFlux'];
paths.fba.fluxAnalysisPath = [paths.fba.fluxPath, filesep, 'fluxAnalysis'];
paths.stats.statPath = [resultPath, filesep, 'resultStatistics'];
% The function initPersephone.m performs the initCobraToolbox.m function, sets the COBRA solver and checks if the required toolboxes are installed (line 1-7). Additionaly it generates the folder structure for the results. IMPORTANT, the output structure used in all 3 tutorials is generated with initPersephone.m. The metadata and the reads table are put through sanity checks. initPersphpne.m ensures that in the metdata the columns with the sample IDs and sample sex have the correct headers and readable data types. It always assumed the first column contains sample IDs. It will create a new updated file and will also update the metadata path accordingly. If you use the test data files you can see that the column sample ID has been changed to ID and the column gender has been altered to Sex. Additionally the data in the column Sex has been changed from M/F to male/female. Important if your own metadata does not have alternative headers accounted for in the function it will raise ahnerror. That is easily fixed by changing your column headers to ID or Sex depending on which data type is causing the issue. Finally we test if the sample IDs in the microbiome data match that of the metadata.
% Create the file structure for the results of all the parts of the
% tutorial (not just this tutorial)
[initialised, statToolboxInstalled, updatedMetadataPath] = initPersephone(resultPath, paths);
% The following ouputs are returned from initPersephone
% initialised - Boolean, indicates if Persephone was succesfully initialised.
% statToolboxInstalled - Boolean, indicates if the statistics toolbox is installed. Parts of Persephone are skipped if false
% updatedMetdataPath  -  The path to the updated metadata file.
% If you look at the directory given in the resultPath variable, you will see the following folder structure.
% -resultDir
%     -HMmodels
%     -personalisedWBMs 
%     -resultFlux
%         -fluxAnalysis
%     -resultMars
%     -resultMgPipe
%     -resultSeqC
%     -resultStatistics
% The content of each of the folders will be discussed in the appropriate sections in the three tutorials.
% We have to also update the metadata path to reflect the updated metadata file.
% paths.General.metadataPath = updatedMetadataPath;
% Now that we have checked that all our software is installed and available, our results directory has been set up and our metadata and read files has been processed and ready for use we can start with processing the metagenomic reads data through MARS.
% 
% Section 2: Running MARS
% In this section we will run MARS to convert the taxonomic reads table to AGORA2 mapped and normalised relative abundance table. MARS works as following:
% First, if the taxonomic assignments and reads tables are separated - it will merge them into a single dataframe
% It removes clade extensions from all taxonomic levels names (e.g. Firmicutes__A & Firmicutes__B will be merged to Firmicutes) to allow for optimal AGORA2 mapping
% It translates certain species based on a pre-defined list to make sure the species names match the ones in the AGORA2 databases
% It removes all reads associated with a taxonomics identification that does not have information up to the specified taxonomic level that is looked at (e.g., looking at species, only reads with information up to the species level will be retained)
% It maps each taxonomic level (kingdom, phylum,species etc.) to the AGORA2 databases. If the taxonomic identification matches it means there is a model present
% The mapped reads are normalised to obtain relative abundances
% All relative abundances under a certain cutoff are removed and the data is re-normalised. The default cutoff is 1e-6.
% Important before running this section is that MARS can also be ran online in your browser on https://mars-pipeline.streamlit.app. The website explains what you have to input and the variables you can use are the same as explained here. Important is that you download the present_species.csv file (file exentsion can differ) and save it in either the 'present' directory in resultMARS named as 'present_species.csv' or give the correct path to MgPipe (discussed in the next section) This ensures that we the MARS output file can still be found by the rest of the functions.
% If we want to run MARS offline we need to make sure all of our dependencies are installed. Follow the code below to install your dependencies and obtain your python path.
% Instructions on how prepare the python environment with 2 options: via miniconda, via user created environment. If you do not have a lot of (python) coding experience the mini-conda way is most user friendly Anaconda: Install Anaconda (https://docs.anaconda.com/miniconda/miniconda-install/)
% open Anaconda Prompt (miniconda3)
%     >> pip install pandas
%     >> pip install numpy
%     >> pip install pyarrow==18.1.0 (latest version of pyarrow gives compatibility issues)
%     >> pip install fastparquet
%     >> pip install openpyxl
% Find the path to the python.exe file in mini-conda(Windows):
% 	-Anaconda Prompt (miniconda3)
% 	-enter "where python" in the prompt
% 	-If no new enivironment is made (the anaconda navigotor only has "base" in environments)
% 	-enter "where python" in the prompt
% 	For macOS and Linux run "which python"
% The python.file exention file location to copy should be in \anaconda3\python.file extension. Paste this path in the pythonPath variable in line 19
% For a user created environment: Make sure you have a working python.exe file on  your computer. Install from https://www.python.org/downloads/. Follow the steps in https://nl.mathworks.com/matlabcentral/answers/1750425-python-virtual-environments-with-matlab. Make sure you install
%     >> pip install pandas
%     >> pip install numpy
%     >> pip install pyarrow
%     >> pip install fastparquet
%     >> pip install openpyxl
% Find the path of the executable of the virtual environment as described
% 		>> import sys
% 		>> sys.executable
% The python.file exention file location to copy should be in \anaconda3\python.file extension. Paste this path in the pythonPath variable in line 19
pythonPath = 'C:\Users\MSPG\miniforge3\python.exe';
% Next we need to clone the MARS repository. Instructions can be found at https://github.com/Unsaif/mars-pipeline. Set the path where the MARS repository was saved as marsRepoPath.
marsRepoPath = 'C:\Users\MSPG\mars-pipeline';

% First we will test if Python is already coupled with Matlab. Then we enter the cloned MARS repository and add the MARS as a readable module. In order to call the apporiate function, we have to enter the MARS directory.
% Check if Python is already coupled to MATLAB, otherwise couple it
pyStatus = pyenv('Version',pythonPath);

% Enter the MARS folder
cd(marsRepoPath);

% Import the entire MARS repository so that the all scripts are on a path
%that MATLAB recognises
MARS = py.importlib.import_module('MARS');

% Enter folder that contains the "main.py" script
cd(strcat(marsRepoPath, '\MARS'));

% Now we can prepare for the actual MARS inputs
% Required inputs
% readsTable - The path to the the taxonomic reads file containing the amount of reads per taxonomic assignment. If the taxonomic assignments are in the reads file the first column needs to be called Taxon. If not the first column of the readsFile and the first column of taxaTable needs to have to same name. See the example files. (we set this already at the start)
% taxaTable - The path to the taxonomic assignments for the taxomomic unit name used in your taxonomic assignment software (e.g., OTU, ASV, OGU). This has to be set to string(missing) if your readsTable already has taxonomic assignment included. Otherwise it need to consists of a column with header similar to the first header in readsTable and the column header "Taxon".
% outputPathMars - The path where the MARS results should be stored. We created this directory in the section 1 and is stored in paths.mars.
% Optional inputs
% sample_read_counts_cutoff - Numeric value for total read counts per sample under which samples are excluded from analysis. Only applies when readsTable contains absolute read counts (not relative abundances). Defaults to 1, with minimum of 1.
% cutoffMars - The cutoff underwich all relative abundances will be considered 0. If smaller relative abundances are kept, numerical difficulties/infeasibilites could occur later on when solving the microbiome models. Also, relative abundances below this value will most likely not have a noticable impact on the flux results as their contribution to fluxes is minimal due to their low relative abundance. Default is 1e-6.
% outputExtensionMars - The file type of your output. Default is csv
% flagLoneSpecies - A boolean, flagLoneSpecies, which is true if the species name does NOT have the genus name already there. False if otherwise. Defaults to false
% taxaSplit - A string,  taxaSplit, which indicates the delimiter used to separate taxonomic levels in the taxonomic assignment.Defaults to '; '
% removeCladeExtensionsFromTaxa - A boolean, removeCladeExtensionsFromTaxa, which removes all reads associated with a taxonomics identification that does not have information up to the specified taxonomic level that is looked at. Defaults to true.
% whichModelDatabase: A string defining if AGORA2, APOLLO, a combination of both or a user-defined database should be used as model database to check presence in.  Allowed Input (case-insensitive): "AGORA2", "APOLLO", "full_db", "user_db". Defaults to "full_db".
% userDatabase_path: A string containing the full path to the user-defined database, which should be in .csv, .txt, .parquet or .xlsx format and have column names = taxonomic levels. Only required if 'whichModelDatabase' is set to "user_db". Note, that the user database needs to have the same structure as the integrated AGORA2 & APOLLO database to function properly!
% We know now which inputs we need to define. Let us start with the required ones.
% Set the path the taxaTable. If you do not have a taxaTable file, put a %
% infront of line 41 and remove % from line 42.
% taxaTable = '';
taxaTable = string(missing)
% The output path stored in the paths variable
outputPathMars = paths.Mars.outputPathMARS; % This is the default path created by initPersephone
% Now let us set the optional inputs. You can change these according to your own dataset, for the test data set we recommend using the settings set here.
% Numeric value for total read counts per sample under which samples are
% excluded from analysis. Only applies when readsTable contains absolute
% read counts (not relative abundances). Defaults to 1, with minimum of 1.
sample_read_counts_cutoff = 1;

% The cutoff value for relative abundances
cutoffMars = 1e-6;

% The file extension for the output files
outputExtensionMars = 'csv';

% The flag if genus name is in the species name
flagLoneSpecies = true;

% The delimeter used to separate taxonomic levels
taxaSplit = ';';

% A boolean specifying if one wants to remove clade name extensions from
% all taxonomic levels of microbiome taxa. If set to false, MARS might find
% significantly less models in AGORA2, as clade extensions are not included
% there.
removeCladeExtensionsFromTaxa = true;

% A string defining if AGORA2, APOLLO, a combination of both or a user-defined
% database should be used as model database to check presence in. 
% Allowed Input (case-insensitive): "AGORA2", "APOLLO", "full_db", "user_db".
% Default: "full_db".
whichModelDatabase="full_db";

% A string containing the full path to the user-defined database,
% which should be in .csv, .txt, .parquet or .xlsx format and
% have column names = taxonomic levels. Only required if whichModelDatabase
% is set to "user_db".
userDatabase_path="";

% We have now defined all our variables for MARS. In order to call the optional arguments, we need to convert them into paired arguments readable by Python through the pyargs function. We will store these converted inputs int he marsOptArg variable.
% Set all optional inputs in Python readable format
marsOptArg = pyargs('cutoff', cutoffMars, 'output_format', string(outputExtensionMars),...
    'flagLoneSpecies',flagLoneSpecies, 'taxaSplit', string(taxaSplit), ...
    'removeCladeExtensionsFromTaxa', removeCladeExtensionsFromTaxa, ...
    'whichModelDatabase', whichModelDatabase, ...
    'sample_read_counts_cutoff', sample_read_counts_cutoff);
% With the optional arguments set we can run MARS
% Run MARS
py.main.process_microbial_abundances(readsTablePath, taxaTable, outputPathMars, marsOptArg);
% return back to the result directory
cd(resultPath);

% Set the relative abundance file path to be used by other functions
relAbunFilePath = [paths.Mars.outputPathMARS, filesep, 'renormalized_mapped_forModelling', filesep, 'renormalized_mapped_forModelling_species.csv'];
% Follow the instructions if you get errors when running line 61:
% ERROR: Unable to resolve the name py.main.process microbial abundances.
%     Check if all you python dependencies are correctly installed. If you had to install them/update them restart MATLAB and rerun the section 1 and section 2 of this tutorial. Check if you have installed pyarrow 18.1.0. To degrade use >> pip uninstall pyarrow. >> pip install pyarrow==18.1.0. It could also be that the MATLAB version is incompatible with the Python version. Either upgrade your MATLAB version or downgrade your Python version.
% ERROR: Python Error: ValueError: Unsupported file type: not found
%     Double check the file extensions of your readsTable and taxaTable inputs. If taxaTable does not exist make sure to set it to string(missing)
% ERROR Python Error: ValueError: Length mismatch: Expected axis has 8 elements, new values have 7 elements or a variation on this
%     Make sure that the taxaSplit is the one that can properly split the taxonomic assignments you have put in. MARS expects that every taxonomical level is defined. If you have only species information in the file, add blank taxonomical information infront in the style of the example (e.g., ; ; ; ; ; ;species name)
% ERROR: Python Error: OSError: Repetition level histogram size mismatch. 
% There is an incompatibility reading parquet files with pyarrow 19.0. This is fixed by degrading to pyarrow 18.1.0. To see the current pyarrow version run
% >> pip show pyarrow
% in the terminal or prompt you used to install Python packages with. If the version is 19.0, degrade by running
% >> pip uninstall pyarrow
% >> pip install pyarrow == 18.1.0

% If MARS finished succesfully then you can find the results in the resultMars directory. 
% metrics - For each taxonomic level, various metric such as alpha and beta diversity are calculated.
% normalized_mapped - For each taxonomic level, only the mapped taxa, abundances are not renormalised and still have the values from the pre-mapping normalisation. Columns do not add up to 1
% normalized_preMapped - For each taxonomic level, all the taxa are normalised to the total reads. Columns add up to 1.
% normalized_unMapped - For each taxonomic level, only the unmapped taxa, abundances are not renormalised and still have the values from the pre-mapping normalisation. Columns do not add up to 1.
% renornmalized_mapped_forModelling - For  each taxonomic level, the mapped taxa are renormalised so that each column adds up to 1.
% The file renormalized_mapped_forModelling_species.csv (file extension can change based on user input) in the ”renormalized_mapped_forModelling” directory is used as input to create community microbiome models. We can expect for about 70-80% of the total reads to be able to be mapped to the AGORA2 database if you use whole-genome shotgun data. Using 16s sequencing, this can drop to 30-40%, mostly because many of the reads do not have information on the species level and are automatically discarded.
% NOTE: If in the normalised_unMappedfile_species file has many species absent it might be worthwile to see if homosynonyms can be found. You can either do this manually or run your list of absent species throught the MARS-ANT pipeline that checks for homosynonyms and presence in the AGORA2 resources. https://mars-pipeline.streamlit.app. Again If the online app is used, please store the files in the appropriate location in your resultDir we defined at the beginning of this part of the tutorial. This allows for seemless integration with the other function and eliminates the need to generate new path variables. If you did fine homosynomys manually, you will need to adjust the names of your taxonomic assignments in either the readsTable or the taxaTable and rerun MARS.
% 
% Section 3: Creating microbiome community models
% In this section we will create community microbiome models using the Microbiome Modelling Toolbox 2.0 [3]. We will use here pan-species models from the AGORA2 resource [2]. Pan-species means that all the strains of a species were combined into one general model that has all the metabolic capabilities of the individual strains. This is done as metagenomic sequencing usually can only distinguish up untill the species level. We will explain how to pipeline generates community metabolic models and how to call on the pipeline (MgPipe). The microbiome modelling toolbox does the following:
% It transforms all microbial reconstructions present in the relative abundance file. The external compartment is changed from [e] to [u] and reactions names are adapted. This is done to prepare the models to exchange with the internal environment of the microbiome [u].
% All metabolites that each microbial reconstruction can exchange with its environment are stored
% All reactions present in each reconstruction are stored as well as their associated subsystems. These will be used to calculate reactionsPresence, reactionAbundance and subsystemAbundance and can be used to describe the metabolic capabilities of indiviual microbioem models.
% All microbe reconstructions that are present in a sample are added together in one big model. The reconstructions can exchange metabolites in the [u] compartment created in the beginning of the function. The microbiome model can exchange metabolites with the outside environemt by transporting metabolites from [u] to [e] and vice versa.
% A microbiome biomass reaction is formulated which looks like: a microbe1 + b micobre2 + c microbe3 + ... The letters a,b and c are the relative abundance for the respective microbes as found in the relative abundance file we calculated via MARS.
% A standard diet is added to the microbiome models and checked if the total growth rate of 1 microbiome exchange per day (read as 1 fecal exchange per day) can be reached.
% Statistics on total amount of reactions, metabolites and microbes for each microbiome are collected and plotted in violin plots.
% To run MgPipe we need to understand the inputs required. NOTE: we here only use a small subset of inputs for MgPipe, if you want a more detailed explanation of the various inputs for the pipeline please look at https://github.com/opencobra/COBRA.tutorials/tree/0d5becf7171844e53d3d43437035917a2bd9cc61/analysis/microbiomeModelingToolbox. Which is also a tutorial in the COBRA toolbox.
% Required inputs
% panModelsPath - The path to the directory with the pan-species models. These can be downloaded from xx
% relAbunFilePath - The path to the file with the AGORA2 mapped relative abundances. This is the present_species.csv file created by MARS in section 2. If the file is stored in the "present" directory in the "resultMars" directory, the variable relAbunFile in the paths variable can be used. Otherwise you have to define it yourself
% computeProfiles - A boolean (true or false) that tells the pipeline if we want to calculate the exctretion and uptake profiles of the microbiome models. We will not do this during this tutorial. However if you are more interested in the microbiome output, you can put this variable to true. This is a time-consuming step and is sped-up by using ibm_cplex as solver. For more information on the extra output if this is set to true look at this tutorial in the COBRA toolbox https://github.com/opencobra/COBRA.tutorials/tree/0d5becf7171844e53d3d43437035917a2bd9cc61/analysis/microbiomeModelingToolbox
% Optional inputs
% numWorkersCreation - Number of workers you want to use for parellelisation. More workers means faster computing times. If you are unsure about the value, follow instructions in the code below. That will help make a decision.
% mgPipeRes - The path where the mgpipe results and the microbiome community models will be stored. The directory is stored in the paths variable created in section 1
% We know how MgPipe works and which inputs we need, we will now define the variables.
% The path to the pan-species directory
panModelsPath = 'D:\OneDrive - National University of Ireland, Galway\panMicrobeModels_Current\Apollo+AGORA2\panSpecies';

% If the relative abundance file is not in the correct MARS folder remove
% the % in front of line 71 and set the correct path
%relAbundFilePath = 'C:Users/User/folder/present/present_species.csv'

% Set computeProfiles variable
computeProfiles = false;

% Set mgPipeRes variable
mgPipeResPath = paths.mgPipe.mgpipePath;

% It can be difficult if you are unfamiliar with the concept to define the number of workers you want to use. Here as some tips you can use. Running "feature('numcores')" will tell you the amount of available cores on your device. An automatic way of chosing your amount of workers is by taking 80% of the available cores on your device which leaves cores available for non-MATLAB related tasks. If working on a cluster or if the device is not required for anythin else, take all available cores. This is recommended when the numCores is greater than 4. If smaller you want want to set the amount of cores to 2 or 3. You cannot go lower than 2 workers. Below is the code that sets the numworkers for you. 
% Find how many available cores you have
numCores = feature('numcores')
% Set the number of workers manually. You can change this value.
numWorkersCreation = 10;

% Remove the % at line 106 if you want to use all cores. This will overwrite
% the workers set in line 104.
% numWorkersCreation = numCores;

% Remove the % at line 108 if you want to use 80% of the cores. This will overwrite
% the workers set in line 104 or 108 if active.
% numWorkersCreation = floor(feature('numCores')*0.8);
% Now that we have set all our inputs we can run MgPipe with the initMgPipe.m function.
% Run MgPipe COPY IN THE COMMAND WINDOW AND RUN FROM THERE. 
% IF RUN IN MLX IT WILL FREEZE YOUR MLX FILE AND YOU WILL NEED TO RESTART MATLAB
% AND COPY AND PASTE THE ENTIRE FILE TO A NEW.MLX VERSION
% initMgPipe(panModelsPath, relAbunFilePath, computeProfiles, "numWorkers",numWorkersCreation, 'resPath', mgPipeResPath, 'solver', solver);
% Warnings might pop up such as
% Individual health status not declared. Analysis will ignore that
% Directory already exists
% The temporary variable 'loadDiet' will be cleared at the beginning of each iteration of the parfor-loop. If 'loadDiet' is used before it is set, a runtime error will occur
% Column headers from the file were modified to make them valid MATLAB identifiers before creating variable names for the table. The original column headers are saved in the VariableDescriptions property. Set 'VariableNamingRule' to 'preserve' to use the original column headers as table variable names.
% Reaction x not in model
% File 'simRes.mat' not found.
% File 'intRes.mat' not found.
% These warnings are expected to occur and can be safely ignored. Press Run Section in the Live Editor to run the MgPipe function.
% If the followig errors occurs, follow the advice below
% ERROR: Modelname cant be found (you dont have the model in the directory or there is some sort of misspelling in the filenames)
% If no errors popped up the microbiome models should now have been created. You can find them in your output folder under mgpipeResults. There you also find the csv files:
% GrowthRates - It gives the maximum growth rate of each microbiome under a rich diet and a pre-defined diet. 
% ModelStatistics - for each microbiome model it gives the total amount of reactions, metabolites and microbes. MicrobiomeModel_sizes.png is a visual representation of this information.
% ModelStatsSummary - gives you the mean, median, min and max of the amount reactions metabolites and microbes considering all microbiome models. 
% ReactionAbundance - gives the relative abundance of a reaction for each metabolic reaction in each microbiome model. Reaction abundances are based on the relative abundance of the microbe and wether or not the microbe has that specific reaction.
% ReactionPresence - says if a reaction is present in at least 1 microbe in a microbiome model.
% SubsystemAbundance - the same as reactionAbundance, however now all reactions are assigned a subsystem and abundances for the same subsystems are summed.
% You will also find .mat files of the microbiome models unconstrained (mirobiota_model_samp_sample ID) or diet constrained (folder Diet, microbiota_model_diet_sample ID). The individual microbe modelsused to create the microbiome models are stored in the modelStorage folder. infeasModels.mat has information on models that could not produce biomass. mapInfo.mat is used to hotstart the pipeline in case there were errors or MATLAB was unexpectly stopped or crashed during execution of the pipeline.

% Section 4: Creation human-microbiome whole-body models (mWBMs)
% In this section we combine the newly created microbiome community models with the previously developed WBMs [4].  We will show the code required to create one mWBM. Later on we will explain which functions to use to generate multiple mWBMs with just one command.
% First we need to load a microbiome model and a WBM model. Here we will load the microbiome model corresponding to sample ID DP_305 and the female WBM.
% Load microbiome model. As we load it into a variable, we also extract the
% actual model
modelM = load([paths.mgPipe.mgpipePath, filesep, 'Diet', filesep, 'microbiota_model_diet_DP305.mat']).model;

% Load the female WBM model 
modelH = loadPSCMfile('Harvetta');

% When combining the WBM with the microbiome model, the [u] section in the microbiome model is change to [luM] which stands for lumen microbiome. Instead of [e] the microbiome now exchanges metabolites with [LuLi] which is the lumen of the large intestine. Reactions are added to ensure metabolites unique to the microbiome model can be taken up from the diet and be excreted in the feces. The function used to combine a microbiome model with a WBM model is called combineHarveyMicrotiota and takes the following inputs
% modelH                        Whole-body metabolic model structure
% modelM                        Microbiome model
% couplingConstraint            Coupling constraint for microbiome model 
% A coupling constraint in this context connects the metabolic activities of the WBM and microbiome model, ensuring that the exchange of metabolites is balanced and physiologically realistic. The value 400 is often used arbitrarily as a scaling factor to adjust the relative metabolic fluxes between the two systems, accounting for their differences in metabolic capacities. This ensures stability and consistency in the model's simulations.
% Some warning might pop up:
%  Reaction EX_biomass[c] is not in model
% The inserted Model contains an old style coupling matrix (A). The Matrix will be converted into a Coupling Matrix (C) and fields will be adapted.
% These are expected and can be safely ignored
% Now we can combine the models
% Set the coupling constraint to 400
couplingConstraint = 400;

% Combine the microbiome mdoel with harvey
mWBM = combineHarveyMicrotiota(modelH, modelM, couplingConstraint);
% Add additional fields to the HM model for more information
mWBM.sex = "female";
mWBM.version = modelH.version;
mWBM.name = modelM.name;
% After the mWBM is created, we set the excretion of miciobiome biomass in the feces to 1 feacal exchange/day [1,1]. This is done to ensure that production of microbiome biomass is the same for all models and that the amount of produced microbiome biomass does not influence the simulation results. We change the constraints of the reaction with the changeRxnBounds function which takes the following inputs:
% model - A metabolic model, in our case mWBM
% rxnNameList - A cell array of reactions that have to have their bounds changed.A cell array of reaction IDs or a string for a single reaction ID that species which reaction(s) need their bounds altered.
% value - A numeric value that say to which value the bound has to be changed to. Can be one value for all reactions, or one value for each reaction to be changed.
% boundType - A string, indicating which bound has to be changed, u= upper, l= lower, b = both.
% We will set the inputs and change the bounds of the faecal microbiome biomass exchange.
% The name of the exchange of faecal microbiome biomass
rxnNameList = 'Excretion_EX_microbiota_LI_biomass[fe]';

% The value of the updated bounds
value = 1;

% Which bounds need to be adjusted, b for both
boundType = 'b';

% Change the bounds of the microbiome biomass fecal excretion to [1,1]
mWBM = changeRxnBounds(mWBM,rxnNameList, value, boundType);

% Here we also define the diet that we want to use to simulate our human-microbiome models with. The default option is EUAverageDiet but other diet options are HighFiberDiet, HighProteinDiet, UnhealthyDiet and VegetarianDiet. If you want to change the diet use one of the afore-mentioned diet names and set that as the Diet variable istead of 'EUAverageDiet'. Using setDietConstraints we put the diet on the mWBM. The function also adds additional metabolites in small amounts (0.1-1 mmol/day). These extra metabolites are added as they are usually not present in the diet files, either pre-made or designed on vmh.life. The metabolites range from ions to small molecules required to ensure microbiome biomass production. The function takes the following inputs
model - model structure, in our case mWBM
diet - Diet option: 'EUAverageDiet' (default)
factor - value between 0 and 1; default is 1, i.e, 100% of the provided diet
% We will define our inputs and set the diet
% Set the chosen diet
diet = 'EUAverageDiet'
factor = 1;
mWBM = setDietConstraints(mWBM, diet, factor);
% Now we have succesfully created a mWBM. The next step is testing if the mWBM is feasible.

% Section 5: Fixing infeasible mWBM
% Now we will check if our create mWBM is feasible. That means can it produce both microbiome biomass and statisfy the body maintanance reaction. First we will set the bound of faecal exchange of microbiome biomass to 1 faecal exchange/day. This was done in line 123, but  sometimes you would already have mWBMs created previously or they are obtained via another person. It is thus good practice to always reset the bounds on the microbiome biomass fecal exchange, even if it seems redudant. 
% The name of the exchange of faecal microbiome biomass
rxnNameList = 'Excretion_EX_microbiota_LI_biomass[fe]';

% The value of the updated bounds
value = 1;

% Which bounds need to be adjusted, b for both
boundType = 'b';

% Change the bounds of the microbiome biomass fecal excretion to [1,1]
mWBM = changeRxnBounds(mWBM,rxnNameList, value, boundType);
% The same goes for the human body maintanance reaction. We will reset it to one to ensure all models are comparable. More information on the 'Whole_body_objective_rxn' can be found in tutorial 2.
mWBM = changeRxnBounds(mWBM, 'Whole_body_objective_rxn', 1, 'b');
% Now that the bounds have been reset we will set the objective of the model to optimise the faecal microbiome biomass exchange. Using changeObjective we can set the reaction we want to have optimised. 
% Set the objective for excretion of microbiome biomass in the feces
mWBM = changeObjective(mWBM,'Excretion_EX_microbiota_LI_biomass[fe]');
% Then we optimise the model with optimizeWBModel. More information on setting objectives and solving models can be found in tutorial 2.
% Solve the HM WBM
solution = optimizeWBModel(mWBM);
% Print the value of the solution
solution.f

% If the sol.f value is 1, great! That means that the model can produce both enough microbiome biomass and can maitain the body maintanance reaction. To ensure we can use the model again we will also save it. The path to where mWBMs are saved was created in section 1 and saved in the paths variable.
% Save the model
save([paths.mWBM.mWBMPath,filesep,'mWBM_DP305_female.mat'],'-struct','mWBM')

% If the output is not 1, it means that the model is unable to create 1 microbiome biomass per day. We will check if the model is able to do it if we open up all dietary reactions. To find all dietary reactions we use the find and contain function to obtain all the reaction indexes that have the prefix Diet_EX_, which indicates dietary exchange reactions. We then use changeRxnBounds to set the lower bound of each dietary exchange reaction to -100000 the unconstrained value used in the WBMs. Then we solve the model again to see if the diet was the problem. The upper bound does not have to be set to 0, as any unused metabolite can pass from the diet to the fecal exchange via various compartments.
% Find all dietary reactions
idx = contains(mWBM.rxns,'Diet_EX');

% Open all dietary exchange reactions
modelHMOpen = changeRxnBounds(mWBM, mWBM.rxns(idx),-1000000,'l');
modelHMOpen.osenseStr = 'max';

% Solve the model
solOpen = optimizeWBModel(modelHMOpen);

% Print the solution value
solOpen.f

% If the solOpen.f value gives 1 that means that means that we can solve the infeasibility by adapting the diet. which is coverd directly below. However if solOpen.f is not 1, the chosen diet is not the problem. You can either a) try another diet or b) check the microbiome explained in the next paragraph.
% More often than not the cause of infeasibility is the absence of micronurtients and can be solved by adding the missing compound to the diet with a value of 0.1 mmol/day. To find out which compound needs to be added we use the function getMissingDietModelHM.m. The function opens up the bounds on all dietary reactions to see if the model can be made feasible through dietary adaptations. If it is solvable with all dietary reactions open, the function will then close random dietary reactions in batches to try and narrow down which reaction is required for feasibility. It will continue to do this untill it has narrowed down a set a essential reactions. As the process is random, running the code twice on the same model might give differing results. If no dietary solution is possible it is good practice to double check that the germ-free (WBM withouth the microbiome) models and the microbiome models are individually feasible with the chosen diet. If the microbiome models are infeasible it might be good to check a) a different diet, b) any irregularities in the present_species file generated from MARS. If a solution cannot be found it is possible to open a issue on the googlegroups for support: https://groups.google.com/g/cobra-toolbox
% The function getMissingDietModelHM takes the following inputs:
% WBM - A WBM model, in our case mWBM
% missingDietComponents - List of diet exchange reactions that are known to be missing or that have been identified in a previous run of this function. We set this to an empty cell array as we have no previous known identified reactions.
% testInitialFeasibility - Boolean to see if an initial feasibility test should be run. We set this to 0 as we already did this beforehand.
% Set inputs
dietExtensions = {};
testInitialFeasibility = 0;
    
% Find missing diet component and add to the model
missingDietComponents = getMissingDietModelHM(mWBM,dietExtensions,testInitialFeasibility);

% Add the missing components to the diet
mWBM.lb(matches(mWBM.rxns,string(missingDietComponents)))=-0.1;

% Save the updated mWBM
save([paths.hmDir,filesep,'mWBM_DP305_female.mat'],'-struct','mWBM');

% Creating multiple mWBMs
% In order to easily generate multiple HM models we need a file that matches the sex with the sample IDs. This will determine if the microbiome model created from that sample will be combined with a female (Harvetta) or male (Harvey) WBM. The diet we will use here is the EUAverageDiet. The metadata path, the path to the microbiome models and the path where the mWBMs should be stored were all defined at the beginning of this tutorial.
% The function createBatchMWBM.m will create for all the microbiome models in the microbiome directory a mWBM. Given that the sample IDs matching to the microbiome models have a corresponding sex in the metadata. The code also checks if the any of the mWBMs are infeasible and then runs the getMissingDietModelHM.m function in a loop to ensure all the mWBMs that are created are feasible and all have the same dietary constraints. This allows for proper comparisons of the flux predictions later on. The function createBatchMWBM.m takes the following inputs:
% microbiomeDir - A string with the path to the directory where the microbiome models are stored, created via MgPipe. Created in section 1
% saveDir - A string with the path to the directory where the mWBMs should be stored. Created in section 1.
% metadataPath - A string with the path to the metadata file. Defined in section 1.
% diet - A string defining which diet should be used. Defaults to EUAverageDiet
% numWorkersCreation - Number of parrelel instances create mWBMs. Same value as used for MgPipe in section 3 is recommended. Defaults to 4.
% numWorkersOptimisation - Number of parrelel instances used to solve the mWBMs, defaults to 2. If the number of workers is too high, it can cause time-out errors. Reducing the number of workers should resolve that issue.
% checkFeasibility -  A boolean to indicate if the function checks if the created mWBMs can grow on the the specified diet if set to true. Defaults to true.
% Now we can set the inputs and run the function.
% Set the diet
diet = 'EUAverageDiet';

% The microbiome path was created already in the beginning of the tutorial
% and stored in the paths variable.
microbiomeDir = paths.mgPipe.mgpipePath;

% The path where the HM models should be stored was already created in the
% beginning of the tutorial and stored in the paths variable.
mWBMdir = paths.mWBM.mWBMPath;

% Set numWorkers if not done so before by removing the % you can alter the
% value
% numWorkersCreation = 1;

% Set checkFeasibility to true
checkFeasibility = true;

% Set the number of workers for optimisation
numWorkersOptimisation = 2;

% Generate multiple HM models
createBatchMWBM(microbiomeDir, mWBMdir, paths.General.metadataPath, "diet",diet, 'numWorkersCreation', numWorkersCreation, 'numWorkersOptimisation', numWorkersOptimisation, "checkFeasibility", checkFeasibility, "solver", solver);
% This concludes this tutorial on how to generate mWBMs. Please see the other tutorials on personalising the WBM models, solving mWBMs models and analysing flux results.
% 
% References (APA style)
% [1] Hulshof, T., Nap, B., Martinelli, F., & Thiele, I. (2024). Microbial abundances retrieved from sequencing data—automated NCBI taxonomy  (MARS): a pipeline to create relative microbial abundance data for the  microbiome modelling toolbox and utilising homosynonyms for efficient  mapping to resources. Bioinformatics Advances, vbae068.
% [2] Heinken, A., Acharya, G., Ravcheev, D. A., Hertel, J., Nyga, M., Okpala, O. E., ... & Thiele, I. (2020). AGORA2: Large scale reconstruction  of the microbiome highlights wide-spread drug-metabolising capacities. BioRxiv, 2020-11.
% [3] Heinken, A., & Thiele, I. (2022). Microbiome Modelling Toolbox 2.0:  efficient, tractable modelling of microbiome communities. Bioinformatics, 38(8), 2367-2368.
% [4] Thiele, I., Sahoo, S., Heinken, A.,  Hertel, J., Heirendt, L., Aurich, M. K., & Fleming, R. M. (2020).  Personalized whole‐body models integrate metabolism, physiology, and the gut microbiome. Molecular systems biology, 16(5), e8982.
% 
