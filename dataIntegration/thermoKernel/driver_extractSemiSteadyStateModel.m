%% input data folder
inputFolder = '~/drive/metaPD/data/metabolomics/inputData';
% specificData
% bibliomicData = 'PD_test1.xlsx';
% bibliomicData = 'PD_test4.xlsx';%allmets
%bibliomicData = 'PD_test5.xlsx';%coremets
bibliomicData = 'PD_test4_withcorerxns.xlsx';%allmets+corerxns
% bibliomicData = 'PD_test5_withcorerxns.xlsx';%coremets+corerxns
% Generic model name
% genericModelName = 'Recon3DModel_301.mat';
genericModelName = 'ReconX.mat';
% results folder
%  resultsFolder = '~/drive/metaPD/results/metabolomics/model/results_weight10000'; %

%% Prepare data for xomics
specificData = preprocessingOmicsModel([inputFolder filesep bibliomicData], 1, 1);

%these weights were too small I think should be minimum -1 for a metabolite or reaction that must be present
specificData.presentMetabolites.weights = specificData.presentMetabolites.weights*100;

%load a decompartmentalised Recon
load('~/drive/metaPD/results/metabolomics/model/modelDecomp_SC.mat')
%%
%load('~/drive/metaPD/results/metabolomics/model/modelSFC.mat')
%% 
metWeights = zeros(size(model.S,1),1);
%omit steady state constraints only for metabolites that are highly
%connected
param.n = 1000; % Connectivity of top 100 metabolites
param.plot = 0; % Do not plot ranked connectivity
param.internal = 1; % Ignore connectivity of stoichiometrically inconsistent part
[rankMetConnectivity,rankMetInd,rankConnectivity] = rankMetabolicConnectivity(model, param);
boolConnected = false(length(metWeights),1);
boolConnected(rankMetInd(1:param.n))=1;
metWeights(boolConnected)=NaN;
% negative weights from clinical data
metWeights = mapAontoB(specificData.presentMetabolites.mets,model.mets,-specificData.presentMetabolites.Entropy,metWeights);
% zero weights for all other metabolites
%metWeights(isnan(metWeights))=0;


%identify the reactions that are exclusively involved in core metabolites, ignoring dummy reactions
coreRxnAbbr = findRxnsFromMets(model,specificData.presentMetabolites.mets);
coreRxnBool = ismember(model.rxns, coreRxnAbbr); 
defaultRxnWeight= ones(size(model.S,2),1)+0.1;
rxnWeights = -defaultRxnWeight;
rxnWeights(~coreRxnBool)=10;

param.printLevel = 1;
param.findThermoConsistentFluxSubset = 0;
param.saveModelSFC = 0;

[thermoModel, thermoModelMetBool, thermoModelRxnBool] = extractSemiSteadyStateModel(model,rxnWeights, metWeights, param);