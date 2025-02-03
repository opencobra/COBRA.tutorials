

load('~/work/sbgCloud/programExperimental/projects/xomics/data/Recon3D_301/Recon3DModel_301_xomics_input.mat')
recon3 = model;
clear model;

load('/home/rfleming/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1/iDopaNeuro1/iDopaNeuro1.mat')
iDN = iDopaNeuro1;
clear iDopaNeuro1;

if ~exist('infeasibleEPproblem','var')
    load('~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/infeasibleEPproblem/infeasibleEPproblem.mat');
end
clear model;

param.printLevel = 1;


[m,nRxns]=size(iDN.S);
n = nnz(iDN.SIntRxnBool);
k = nnz(~iDN.SIntRxnBool);

% case 'pdco'
%     %constraint matrix
%     EPproblem.A  =[...
%         N,     -N,    Omn,     B;
%         In,    -In,    -In,   Onk;
%         C,     -C,    Ocn,     D];


param.excludedReactionLB = false(size(EPproblem.A,2),1);
param.excludedReactionLB(1:2*n) = 1;

param.excludedReactions = false(size(EPproblem.A,2),1);
param.excludedReactions(2*n+1:3*n+k) = 1;
%param.excludedReactions(m+1:end) = 1;

param.steadyStateRelax = 0;

[solution, relaxedEPproblem] = relaxedFBA(EPproblem, param);

%%
%lower bounds on reactions
boolRelaxLower = solution.p~=0;

find(boolRelaxLower)
EPproblem.lb(boolRelaxLower)
relaxedEPproblem.lb(boolRelaxLower)



if any(boolRelaxLower)
    
    boolRxnRelaxLowerFwd = false(size(iDN.S,2),1);
    boolRxnRelaxLowerFwd(iDN.SIntRxnBool) = boolRelaxLower(1:n,1);
    
    printConstraints(iDN,-inf,inf,boolRxnRelaxLowerFwd)
    
    lowerTfwd = table(...
        iDN.rxns(boolRxnRelaxLowerFwd),...
        recon3.lb(ismember(recon3.rxns,iDN.rxns(boolRxnRelaxLowerFwd))),...
        iDN.lb(boolRxnRelaxLowerFwd),...
        EPproblem.lb(boolRelaxLower(1:n,1)),...
        relaxedEPproblem.lb(boolRelaxLower(1:n,1)),...
        'VariableNames',{'rxns','recon','iDN','EP_lb','relaxEP_lb'})
    
    boolRxnRelaxLowerRev = false(size(iDN.S,2),1);
    boolRxnRelaxLowerRev(iDN.SIntRxnBool) = boolRelaxLower(n+1:2*n,1);
    
    printConstraints(iDN,-inf,inf,boolRxnRelaxLowerRev)
    
    lowerTrev = table(...
        iDN.rxns(boolRxnRelaxLowerRev),...
        recon3.lb(ismember(recon3.rxns,iDN.rxns(boolRxnRelaxLowerRev))),...
        iDN.lb(boolRxnRelaxLowerRev),...
        EPproblem.lb(boolRelaxLower(1:n,1)),...
        relaxedEPproblem.lb(boolRelaxLower(1:n,1)),...
        'VariableNames',{'rxns','recon','iDN','EP_lb','relaxEP_lb'})
end

%%
%upper bounds on reactions
boolRelaxUpper = solution.q~=0;

find(boolRelaxUpper)
EPproblem.ub(boolRelaxUpper)
relaxedEPproblem.ub(boolRelaxUpper)

if any(boolRelaxUpper)
    
    boolRxnRelaxUpperFwd = false(size(iDN.S,2),1);
    boolRxnRelaxUpperFwd(iDN.SIntRxnBool) = boolRelaxUpper(1:n,1);
    
    printConstraints(iDN,-inf,inf,boolRxnRelaxUpperFwd)
    
    upperTfwd = table(...
        iDN.rxns(boolRxnRelaxUpperFwd),...
        recon3.ub(ismember(recon3.rxns,iDN.rxns(boolRxnRelaxUpperFwd))),...
        iDN.ub(boolRxnRelaxUpperFwd),...
        EPproblem.ub(boolRelaxUpper(1:n,1)),...
        relaxedEPproblem.ub(boolRelaxUpper(1:n,1)),...
        'VariableNames',{'rxns','recon','iDN','EP_ub','relaxEP_ub'})
    
    boolRxnRelaxUpperRev = false(size(iDN.S,2),1);
    boolRxnRelaxUpperRev(iDN.SIntRxnBool) = boolRelaxUpper(n+1:2*n,1);
    
    printConstraints(iDN,-inf,inf,boolRxnRelaxUpperRev)
    
    upperTrev = table(...
        iDN.rxns(boolRxnRelaxUpperRev),...
        recon3.ub(ismember(recon3.rxns,iDN.rxns(boolRxnRelaxUpperRev))),...
        iDN.ub(boolRxnRelaxUpperRev),...
        EPproblem.ub(boolRelaxUpper(n+1:2*n,1)),...
        relaxedEPproblem.ub(boolRelaxUpper(n+1:2*n,1)),...
        'VariableNames',{'rxns','recon','iDN','EP_ub','relaxEP_ub'})
end


printConstraints(iDN,-inf,inf,iDN.lb<0 & iDN.ub==0 & iDN.SIntRxnBool)





