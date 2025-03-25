% Bioinformatic processing of sequencing data
% Authors: Wiley Barton
% Created: 2025.02.04
% Requirements:
% To successfully follow this tutorial, you need to have the following installed on your system:
%     Docker Desktop or Engine (tested@ 4.37.1 (178610))
%     MATLAB (tested@R2024b)
%    SeqC Repository (Bioinformatics pipeline for taxonomic classification)
% Introduction
% This tutorial is part of a series to support the use of the Persephone pipeline. This tutorial goes through the steps of the overall pipeline that interact directly with the Sequence Conversion (SeqC) pipeline. SeqC is a Docker-based bioinformatic environment with the purpose of providing a standardised, efficient, and portable means to generate the microbial taxonomic inputs for the rest of Persephone from raw sequencing data. SeqC performs quality control of fastq files with Kneaddata [1], taxonomic assignment to reads with a combination of Kraken2 [2] and bracken [3]--using custom assignment databases derived from AGORA2 [4] and or APOLLO [5], and additional configuration with MARS [6].
% Section 1: Environment Preparation
% The SeqC pipeline is executed using the runSeqC function, which is called by the main runPersephone function. This function relies on configurations set in the configPersephone.m file.
% Downloading SeqC
% SeqC is included in the COBRA Toolbox. If you have not installed it yet, refer to previous tutorials or follow the instructions available here: https://github.com/opencobra/cobratoolbox
% Alternatively, you can manually clone the SeqC repository by running the following command from you systems command line:
%git clone git@gitlab.com:thielelab/wbm_modelingcode.git --branch master

% Configuring SeqC Paths
% Once SeqC is installed, you need to specify key file paths in MATLAB:
% Character array variable specifying the folder where seqC repository is stored
% e.g., 'C:\Users\cobratoolbox\src\analysis\persephone\SeqC_pipeline'
paths.seqC.repoPathSeqC = '';
% Character array variable specifying the folder
% where the final output of SeqC is stored.
resultPath = '';


paths.seqC.outputPathSeqC = [resultPath filesep, 'ResultSeqC'];
% *REQUIRED*. Character array variable of the file name
% containing sample IDs for FASTQ files (e.g., sample_id.txt)
paths.seqC.fileIDSeqC = 'sample_id_demo.txt';
% These variables control how SeqC processes sequencing data and determines available computational resources.

% Set the paths to all inputs for seqC
% Logical variable indicating if intermediary outputs are
% retained (e.g., post-QC FASTQ). False results in the
% singular output of MARS, and the deletion of all
% intermediary content once SeqC completes.
paths.seqC.procKeepSeqC = false;
% Numeric, the maximum amount of memory allowed in gigabytes.
paths.seqC.maxMemSeqC = 20;

% Numeric, the maximum number of threads allowed.
paths.seqC.maxCpuSeqC = 4;

% Numeric, the maximum number of processes allowed.
paths.seqC.maxProcSeqC = 4;
% Logical variable indicating if additional debugging
% messaging should be included in the log.
paths.seqC.debugSeqC = false;

% Set the paths to all inputs for seqC
%% 1.5 MARS inputs
% Logical variable indicating if taxonomic mapping by MARS
% is performed independently of SeqC.
paths.Mars.flagMars = true;
% *REQUIRED*. Character array variable with path to the microbiome taxonomy
% and read abundance file.
paths.Mars.readsTablePath = '';
% Character array variable to the folder where the output of MARS is stored.
paths.Mars.outputPathMars = [resultPath filesep, 'ResultMars'];

% Character array variable indicating the desired file format for saving
% outputs: this is needed for the relative abundance file path used in the
% next line
paths.Mars.outputExtensionMars = 'csv';

% Path to relative abundance file
paths.Mars.relAbunFilePath = [paths.Mars.outputExtensionMars filesep 'present' filesep 'present_species' filesep paths.Mars.outputExtensionMars];

% Numeric value for total read counts per sample under which samples are
% excluded from analysis. Only applies when readsTable contains absolute
% read counts (not relative abundances). Defaults to 1, with minimum of 1.
paths.Mars.sample_read_counts_cutoff = 1;

% Numeric value under which relative abundances are considered to be zero
% default = 1e-6
paths.Mars.cutoffMars = 1e-6;
% String to the file where OTUs are matched to taxonomic assignments.
% OPTIONAL if the taxonomic assignments are already in the readsTable.
% REQUIRED if not.
paths.Mars.OTUTable = "";

% A boolean to indicate if the genus name is in the name of the species e.g.
% Prevotella copri. If genus name is in species name set to false.
% Otherwise set to true. OPTIONAL, defaults to false.
paths.Mars.flagLoneSpecies = true;

% The delimeter used to separate taxonomic levels
paths.Mars.taxaDelimiter = ';';

% A boolean specifying if one wants to remove clade name extensions from
% all taxonomic levels of microbiome taxa. If set to false, MARS might find
% significantly less models in AGORA2, as clade extensions are not included
% there.
paths.Mars.removeClade = true;
% A string defining if AGORA2, APOLLO, a combination of both or a user-defined
% database should be used as model database to check presence in.
% Allowed Input (case-insensitive): "AGORA2", "APOLLO", "full_db", "user_db".
% Default: "full_db".
paths.Mars.reconstructionDb = "full_db";
% A string containing the full path to the user-defined database,
% which should be in .csv, .txt, .parquet or .xlsx format and
% have column names = taxonomic levels. Only required if reconstructionDb
% is set to "user_db".
paths.Mars.userDbPath  = "";
paths.Mars.taxaTable = "";

% Section 2: Running SeqC
% The function runSeqC begins by setting the working directory to that of the SeqC repository, and assessing environmental variables and ensuring they are correctly formatted.
% Now that the environment is set up, you can execute SeqC by calling:
runSeqC(...
            paths.seqC.repoPathSeqC, ...
            paths.seqC.outputPathSeqC, ...
            paths.seqC.fileIDSeqC, ...
            paths.seqC.procKeepSeqC, ...
            paths.seqC.maxMemSeqC, ...
            paths.seqC.maxCpuSeqC, ...
            paths.seqC.maxProcSeqC, ...
            paths.seqC.debugSeqC, ...
            paths.Mars.readsTablePath, ...
            paths.Mars.outputPathMars, ...
            paths.Mars.outputExtensionMars, ...
            paths.Mars.relAbunFilePath, ...
            paths.Mars.sample_read_counts_cutoff, ...
            paths.Mars.cutoffMars, ...
            paths.Mars.OTUTable, ...
            paths.Mars.flagLoneSpecies, ...
            paths.Mars.taxaDelimiter, ...
            paths.Mars.removeClade, ... 
            paths.Mars.reconstructionDb, ...
            paths.Mars.userDbPath, ... 
            paths.Mars.taxaTable ...
    );

% Checking Your Operating System
% The tutorial automatically detects your OS and ensures that Docker is properly configured:
%% Determine Operating System
% Some arguments are OS specific
if ismac
    vOS = 'mac';
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']); % Ensure Docker is found
elseif isunix
    vOS = 'unix';
elseif ispc
    vOS = 'win';
else
    error('Unsupported operating system.');
end
% Ensure Docker is running before executing SeqC, especially on Windows.

% Estimating Disk Usage
% SeqC requires a considerable amount of storage. To estimate space requirements, an optional estimation of the size of outputs generated is calculated below.
%% Determine directory size and estimate usage exapansion
dirPath = fullfile(paths.seqC.repoPathSeqC,'seqc_input/'); % Directory path
totalBytes = getDirectorySize(dirPath); % Function from previous response
totalMB = totalBytes / (1024^2); % Convert to MB
totalGB = totalBytes / (1024^3); % Convert to GB
inflateRatio = 3.2; % inflation term
inflateGB = totalGB * inflateRatio;
msgDsize = sprintf('Total size of directory: %.2f GB\nExpected inflation size: %.2f GB\n', totalGB, inflateGB);


% Before running SeqC, we need to define and initialize key directories. This ensures that the pipeline can correctly locate input files and store output results.
% In MATLAB, initialize the required paths with:
%% Initialize Paths
vdir_init = cd;
vdir_out_seqc = 'seqc_output';
vdir_out_mars = fullfile(vdir_out_seqc, 'mars_out');
% Set system to seqc repo
cd(paths.seqC.repoPathSeqC);

% Before passing computational resource constraints to the Docker container, we convert numeric values to strings for compatibility. 
%% Convert Numeric Inputs to Strings
maxCpuSeqC = num2str(paths.seqC.maxCpuSeqC);
maxMemSeqC = num2str(paths.seqC.maxMemSeqC);
maxProcSeqC = num2str(paths.seqC.maxProcSeqC);
if isnumeric(paths.Mars.cutoffMars)
    paths.Mars.cutoffMars = sprintf('%f', paths.Mars.cutoffMars);
end
% Then, following command constructs Docker build arguments using the previously converted string values.
%% Build Docker Options
% Hardware params
comm_build_opt_hw = sprintf('--build-arg varg_cpu_max=%s --build-arg varg_mem_max=%s --build-arg varg_proc_max=%s', ...
                            paths.seqC.maxCpuSeqC, paths.seqC.maxMemSeqC, paths.seqC.maxProcSeqC);

% To ensure MARS settings are correctly passed to Docker, we construct a command string with the necessary arguments. This configures MARS within the SeqC pipeline for proper execution inside the container.
% MARS params
comm_build_opt_mars = sprintf(['--build-arg varg_mars_outputExtensionMars=%s' ...
' --build-arg varg_mars_sample_read_counts_cutoff=%d' ...
' --build-arg varg_mars_cutoffMars=%s' ...
' --build-arg varg_mars_flagLoneSpecies=%s' ...
' --build-arg varg_mars_taxaDelimiter="%s"' ...
' --build-arg varg_mars_removeClade =%s' ...
' --build-arg varg_mars_reconstructionDb=%s'], ...
paths.Mars.outputExtensionMars, paths.Mars.sample_read_counts_cutoff, paths.Mars.cutoffMars, ...
string(paths.Mars.flagLoneSpecies), string(paths.Mars.taxaDelimiter), string(paths.Mars.removeClade), ...
paths.Mars.reconstructionDb);

% Optional MARS parameters are appended to the Docker build command only if they are provided. This ensures a flexible and efficient configuration within the SeqC pipeline.
%% Append Optional Build Arguments - MARS
% exclude if empty
optionalParams = {paths.Mars.OTUTable, paths.Mars.readsTablePath, paths.Mars.relAbunFilePath, paths.Mars.userDbPath, paths.Mars.taxaTable};
paramNames = {'varg_mars_OTUTable', 'varg_mars_readsTablePath', 'varg_mars_relAbunFilePath', 'varg_mars_userDbPath ', 'varg_mars_taxaTable'};
for vi = 1:length(optionalParams)
    if ~ismissing(optionalParams{vi})
        if ~isempty(optionalParams{vi}) && ~strcmpi(optionalParams{vi}, "")
            comm_build_opt_mars = sprintf('%s --build-arg %s=%s', comm_build_opt_mars, paramNames{vi}, optionalParams{vi});
        end
    end
end


% The Docker image is built with hardware and MARS-specific options, ensuring proper resource allocation. The run command is configured for interactive and non-interactive execution, adapting to system constraints.
%% Build Docker Image command
comm_build = sprintf('docker build -t dock_seqc --ulimit nofile=65536:65536 %s %s .', comm_build_opt_hw, comm_build_opt_mars);

%% Docker run commands
% core run command
comm_run_core = 'docker run --interactive --tty --user 0 --rm --mount';
% sans interactive
comm_run_core = sprintf('docker run --tty --user 0 --rm --memory=%s --cpus=%s --mount',sprintf('%sg',maxMemSeqC),maxCpuSeqC);


% The user can select a taxonomic database (AGORA, APOLLO, or combined) based on user input and constructs a command to run the database creation script with the chosen database and a human contamination filter.
%% Set Database Assignment Command
switch paths.Mars.reconstructionDb
    case 'AGORA'
        comm_run_db_kb = '-s "tool_k2_agora"';
    case 'APOLLO'
        comm_run_db_kb = '-s "tool_k2_apollo"';
    case 'full_db'
        comm_run_db_kb = '-s "tool_k2_agora2apollo"';
    otherwise
        comm_run_db_kb = '-s "tool_k2_agora2apollo"'; % Default case
end
comm_run_db_kd = '-s "host_kd_hsapcontam"';
comm_run_db = sprintf('BASH_seqc_makedb.sh %s %s', comm_run_db_kd, comm_run_db_kb);


% Next, we construct the command to run the SeqC pipeline, optionally appending flags for debugging and keeping intermediate files based on user settings. We then format the full command with input data directory, sample IDs, and other required parameters.
%% Construct Command for Running SeqC
comm_mama_help = 'BASH_seqc_mama.sh -h';
comm_mama_full = 'BASH_seqc_mama.sh';
% append optional flags
if paths.seqC.debugSeqC
    comm_mama_full = [comm_mama_full ' -b'];
end
if paths.seqC.procKeepSeqC
    comm_mama_full = [comm_mama_full ' -k'];
end
comm_mama_full = sprintf('%s -i "step0_data_in/" -n "%s" -r "SR" -s 0', comm_mama_full, paths.seqC.fileIDSeqC);

% Next, we construct the command to run the SeqC process in a Docker container, adjusting the volume and directory mappings based on the operating system (Unix, Mac, or Windows). This binds input and output directories and specifies where the processing data will be stored within the container.
% Append volume mapping commands to core
% OS sensitive
if strcmp(vOS, 'unix')
    comm_run_main = sprintf('%s "type=bind,src=$(pwd)/seqc_input,target=/home/seqc_user/seqc_project/step0_data_in" --mount "type=bind,src=$(pwd)/seqc_output,target=/home/seqc_user/seqc_project/final_reports" --mount "type=volume,dst=/DB,volume-driver=local,volume-opt=type=none,volume-opt=o=bind,volume-opt=device=$(pwd)/seqc_proc" dock_seqc /bin/bash', comm_run_core);
%    comm_exit_mv = 'mv -r $(pwd)/seqc_proc/DEPO_proc/* $(pwd)/seqc_output'
elseif strcmp(vOS, 'mac')
    comm_run_main = sprintf('%s "type=bind,src=$(pwd)/seqc_input,target=/home/seqc_user/seqc_project/step0_data_in" --mount "type=bind,src=$(pwd)/seqc_output,target=/home/seqc_user/seqc_project/final_reports" --mount "type=volume,dst=/DB,volume-driver=local,volume-opt=type=none,volume-opt=o=bind,volume-opt=device=$(pwd)/seqc_proc" dock_seqc /bin/bash', comm_run_core);
%    comm_exit_mv = 'mv -r $(pwd)/seqc_proc/DEPO_proc/* $(pwd)/seqc_output'
elseif strcmp(vOS, 'win')
    comm_run_main = sprintf('%s "type=bind,src=%s\\seqc_input,target=/home/seqc_user/seqc_project/step0_data_in" --mount "type=bind,src=%s\\seqc_output,target=/home/seqc_user/seqc_project/final_reports" --mount "type=bind,src=%s\\seqc_proc,target=/DB" dock_seqc /bin/bash', comm_run_core, pwd, pwd, pwd);
%    comm_exit_mv = 'mv -r .\seqc_proc\DEPO_proc\* .\seqc_output\'
end

% Section 3: Running Docker for SeqC
% SeqC runs within a Docker container. The following commands:
% Build the Docker image if it does not exist
% Run the SeqC pipeline
% Step 1: Build the Docker Image
% Once all the variables and Docker statements are constructed, Docker is engaged and an image of SeqC is created if a previous image is not found.
% check for preexisting image
    imageName = 'dock_seqc';
[status, cmdout] = system(['docker images -q ' imageName]);

if isempty(strtrim(cmdout))
    disp(['Image "' imageName '" does NOT exist. Now creating...']);
    disp(' > Building SeqC docker image, wait time ~10min.');
    [status, cmdout] = system(comm_build);
    if status ~= 0, error('Docker build failed:\n%s', cmdout); end
else
    disp(['Docker Image "' imageName '" exists.']);
end

% Step 2: Run the Pipeline
% Once the Docker image is built, the following commands that will:
% Test to confirm image viability
% Establish required databases
% Proccess sequencing files
% First we run a test of the MAMA script to verify if the Docker image and related commands work correctly. If it fails, a warning is displayed.
% Test MAMA script
[status, cmdout] = system(sprintf('%s %s',comm_run_main, comm_mama_help));
if status ~= 0, warning('MAMA test failed:\n%s', cmdout); end

% Next we must initiate the database setup for the pipeline, with an estimated wait time of 30 minutes.
% Run database creation
disp(' > Running database setup, wait time ~30min....');
[status, cmdout] = system(sprintf('%s %s',comm_run_main, comm_run_db));
if status ~= 0, error('Database setup failed:\n%s', cmdout); end

% And finaly, this command starts the full SeqC processing pipeline:
% Run full SeqC pipeline
disp(sprintf(' > SeqC Processing Begins...\n%s',msgDsize));
[status, cmdout] = system(sprintf('%s %s',comm_run_main, comm_mama_full));
if status ~= 0, error('SeqC pipeline execution failed:\n%s', cmdout); end
disp(' > SeqC Processing Ends.');

% Section 4: Managing Output Files
% After processing, results are stored in the outputPathSeqC folder. If MARS is used, additional outputs are placed in outputPathMars, this file structure is important for advancement in the Persephone workflow
% Move final output
movefile(fullfile(vdir_out_seqc, '*'), paths.seqC.outputPathSeqC);
% Update mars path
vdir_out_mars = fullfile(paths.seqC.outputPathSeqC, 'mars_out');
if ~strcmp(paths.seqC.outputPathSeqC, paths.Mars.outputPathMars)
    movefile(fullfile(vdir_out_mars, '*'), paths.Mars.outputPathMars);
    vdir_out_mars = fullfile(paths.seqC.outputPathSeqC);
end
