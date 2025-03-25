initCobraToolbox();

% Can also be different versions
changeCobraSolver('ibm_cplex');


pe = pyenv;
if pe.Version == ""
    error("Python not installed, required for running MARS")
end
% Do we need a check for Pandas, JSON or Numpy? - required for MARS

% Are all paths set up correctly

microbiomeModels = '';
readsFilePath = '';
OTUFilePath = '';

paths2Check = {microbiomeModels; dietFilePath; readsFilePath; OTUFilePath};

for i = 1:size(paths2Check,1)
    if i ==1 
        if ~isfolder(paths2Check{i})
            error(concat(microbiomeModels{i}, ' does not seems to exist on your system. Please check if your path is formulated correctly'))
        end
    else
        if ~isfile(paths2Check{i}) && ~isempty(paths2Check{i})
            error(concat(microbiomeModels{i}, ' does not seems to exist on your system. Please check if your path is formulated correctly'))
        end
    end
end
