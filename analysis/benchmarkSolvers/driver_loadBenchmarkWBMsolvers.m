%% Choose between different models. 
% Add your own via the switch statement below.

if ~exist(modelToUse,'var')
    %modelToUse  = 'iCoreED'; % medium - TODO could not compute extremePools -find out why.
    %modelToUse = 'iDopaNeuroC';
    % modelToUse = 'Recon3T';
    modelToUse = 'Harvey';
    %modelToUse = 'Harvetta';
end

%% Set local parameters
% If true, it will relax tight bounds
relaxTightBounds=0;
lowerExponent = 4; %the minimum difference between ub_j and lb_j is 10^(lowerExponent)
higherExponent = 10;

% Load the selected model and make any adjustments necessary.
switch modelToUse
    case 'iCoreED'
        if ~exist('iCoreED','var')
            load ~/drive/sbgCloud/projects/variationalKinetics/data/iCoreED/iCoreED_modelT.mat
        end
        model    = modelT;

        if ~isfield(model,'biomassRxnAbbr')
            model.biomassRxnAbbr='Biomass_Ecoli_core_w_GAM';
        end
    case 'iDopaNeuroC'
        load('~/drive/sbgCloud/projects/variationalKinetics/data/iDopaNeuro/iDopaNeuroC.mat')
        %load('~/drive/sbgCloud/data/models/published/iDopaNeuro/iDopaNeuroC.mat')
        model = iDopaNeuroC;
        model.description='iDopaNeuroC';
        model.version='1.0';
    case 'Recon3T'
        load('~/drive/sbgCloud/projects/variationalKinetics/data/Recon3T/Recon3DModel_301_xomics_input.mat');
    case 'Harvey'
        if isfile('~/drive/sbgCloud/code/wbm_modelingcode/WBM_reconstructions/Harvey_1_04c_lifted.mat')
            load('~/drive/sbgCloud/code/wbm_modelingcode/WBM_reconstructions/Harvey_1_04c_lifted.mat')
            model=male;
            model.osenseStr='max';
        else
            load('~/drive/sbgCloud/code/wbm_modelingcode/WBM_reconstructions/Harvey_1_04c.mat')
            %load('~/drive/sbgCloud/projects/variationalKinetics/data/WBM/Harvey_1_04c.mat')
            %load('~/drive/sbgCloud/data/models/published/Harvey_Harvetta/20191104_Harvey_1_01c.mat')
            model=male;
            model.osenseStr='max';
            model = changeObjective(model,model.rxns(contains(model.rxns,'Whole')));

            model = homogeniseCouplingConstraints(model);

            %% Reformulate coupling constraints
            % Using hierarchical lifting as described here
            % https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-240
            [male] = liftCouplingConstraints(model);
            save('~/drive/sbgCloud/code/wbm_modelingcode/WBM_reconstructions/Harvey_1_04c_lifted.mat','male');
            return
        end

    case 'Harvetta'
        if isfile('~/drive/sbgCloud/code/wbm_modelingcode/WBM_reconstructions/Harvetta_1_04c_lifted.mat')
            load('~/drive/sbgCloud/code/wbm_modelingcode/WBM_reconstructions/Harvetta_1_04c_lifted.mat')
            model=female;
            model.osenseStr='max';
        else
            load('~/drive/sbgCloud/code/wbm_modelingcode/WBM_reconstructions/Harvetta_1_04c.mat')
            %load('~/drive/sbgCloud/projects/variationalKinetics/data/WBM/Harvetta_1_04c.mat')
            %load('~/drive/sbgCloud/data/models/published/Harvey_Harvetta/20191104_Harvetta_1_01c.mat')
            model=female;
            model.osenseStr='max';
            model = changeObjective(model,model.rxns(contains(model.rxns,'Whole')));

            model = homogeniseCouplingConstraints(model);
            
            %% Reformulate coupling constraints
            % Using hierarchical lifting as described here
            % https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-240
            [female] = liftCouplingConstraints(model);
            save('~/drive/sbgCloud/code/wbm_modelingcode/WBM_reconstructions/Harvetta_1_04c_lifted.mat','female');
            return
        end
    otherwise
        errror('Unrecognised modelToUse')
end



