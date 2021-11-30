
tutorialPath = fileparts(which('tutorial_template.mlx'));


fileList = dir(fullfile(tutorialPath, ['**' filesep '*.*']));  %get list of files and folders in any subfolder

fileList = struct2table(fileList);

mlxFileBool = contains(fileList.name,'.mlx');


fileList = fileList.name(mlxFileBool)