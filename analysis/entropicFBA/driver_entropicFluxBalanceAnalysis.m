%load a model and compute a solution to the bounded entropy problem using

param.internalBounds='';
param.internalBounds='directional';
param.feasTol=1e-7;

param.method='fluxesConcentrations';

%param.solver='pdcoPrimal';
%param.solver='mskexpcone';
%solver='mskenoptPrimal';
param.solver='pdco';
%param.solver='mosek';


if strcmp(param.solver,'mosek')
    %set default mosek parameters for this type of problem
    param=mosekParamSetEFBA(param);
end

param.printLevel=1;

basePath='~/work/sbgCloud';

resultsDirectory=[basePath '/programModelling/projects/thermoModel/results/entropicFBA'];

if ~exist('model','var')
    modelToUse='iDN';
    %modelToUse='ecoli_core';
    %modelToUse = 'fork';
    %modelToUse = 'single';
    
    switch modelToUse
        case 'single'
            driver_singleModel
        case 'fork'
            driver_forkModel
            
        case 'ecoli_core'
            dataDirectory=[basePath '/data'];
            load([dataDirectory '/models/published/e_coli_core.mat']);
            
        case 'recon3'
            genericModelName = 'Recon3DModel_301_thermoFeas_BBF.mat';
            inputFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics'];
            load([inputFolder filesep genericModelName])
            %load([dataDirectory '/models/unpublished/2016_07_13_Recon3_SConsistencyCheck.mat']);
            %model=modelRecon3model;

            
            %change onto uMol/gDW/hr
            model.lb = model.lb;
            model.ub = model.ub;
        case 'iDN'
            load('~/work/sbgCloud/programExperimental/projects/xomics/results/iDN1/gurobi_inputData_oneRxnsPerActiveGene_transcriptomicsT0_mediaBBRelaxed_closedIons_thermoKernel_noBPL_inactiveGenesT_limit10k/Model.mat')
            model = Model;
            clear Model
            clear modelGenerationReport
    end
end

%projection
%model.S=N;
%[PR,PN,PC,PL]=subspaceProjector(model,printLevel-1,'all',1);


%reaction parameters
cf=0;
cr=0;
g='two';

%metabolite parameters
u0='zero';
f=1;

param.maxElementaryFlux=inf;
param.internalBounds='original';
    
[solution, model] = entropicFluxBalanceAnalysis(model,'min',cf,cr,g,u0,f,param);

[v,vf,vr,vt,lambda,zeta,eta,alpha,beta,stat]=...
    deal(solution.v,solution.vf,solution.vr,solution.vt,solution.lambda,solution.zeta,solution.eta,solution.alpha,solution.beta,solution.stat);
    
g=deal(model.g);
if isfield(solution,'x')
    [x, x0, delta, epsilon, eta] = deal(solution.x, solution.x0, solution.delta, solution.epsilon, solution.eta);
    [u0,f]=deal(model.u0,model.f);
    dx = x - x0;
end
    
if exist('osenseStr', 'var') 
    if isempty(osenseStr)
        model.osenseStr = 'max';
    else
        model.osenseStr = osenseStr;
    end
else
    if isfield(model, 'osenseStr')
        model.osenseStr = model.osenseStr;
    else
        model.osenseStr = 'max';
    end
end
[~,osense] = getObjectiveSense(model);


N = model.S(:,model.SConsistentRxnBool);
B = model.S(:,~model.SConsistentRxnBool);
if isfield(model,'C')
    C = model.C(:,model.SConsistentRxnBool);
    nabla = solution.nabla;
end

e   = ones(length(vf),1);
ci = osense*model.c(model.SConsistentRxnBool);

if 0
   %get nullspace of N
    [Z,rankS]=getNullSpace(N,0);
else
    Z=[];
end
%%
if param.printLevel>2 | 1
    detailed=1;
else
    detailed=0;
end

fprintf('\n%s\n',' ------ driver_entropicFluxBalanceAnalysis ---------')
fprintf('\n%s\n','Biochemical optimality conditions')
if isfield(solution,'x')
    fprintf('%8.2g %s\n',norm(N*(vf - vr) - x + x0 - model.b,inf), '|| N*(vf - vr) - x + x0 - b ||_inf');
else
    fprintf('%8.2g %s\n',norm(model.S*v - model.b,inf), '|| S*v - b ||_inf');
end
if detailed
    fprintf('%8.2g %s\n',sum(vf)+sum(vr),'sum(vf)+sum(vr) >= 0');
    fprintf('%8.2g %s\n',norm(g.*reallog(vf),inf),'|| g*log(vf) ||_inf');
    fprintf('%8.2g %s\n',norm(g.*reallog(vr),inf),'|| g*log(vr) ||_inf');
    fprintf('%8.2g %s\n',norm(cr,inf),'|| cr ||_inf');
    fprintf('%8.2g %s\n',norm(cf,inf),'|| cf ||_inf');
    fprintf('%8.2g %s\n',norm(ci,inf),'|| ci ||_inf');
    fprintf('%8.2g %s\n',norm(N'*lambda,inf),'|| N''*lambda ||_inf');
    fprintf('%8.2g %s\n',norm(zeta,inf),'|| zeta ||_inf');
    fprintf('%8.2g %s\n',norm(alpha,inf),'|| alpha ||_inf');
    fprintf('%8.2g %s\n',norm(beta,inf),'|| beta ||_inf');
    if isfield(model,'C')
        fprintf('%8.2g %s\n',norm(C'*nabla,inf),'|| C''*nabla ||_inf');
    end
end

if isfield(model,'C')
    fprintf('%8.2g %s\n',norm(g.*reallog(vf) + ci + cf + N'*lambda + C'*nabla + zeta + alpha,inf), '|| g*log(vf) + ci + cf + N''*lambda + C''*nabla + zeta + alpha ||_inf');
    fprintf('%8.2g %s\n',norm(g.*reallog(vr) - ci + cr - N'*lambda - C'*nabla - zeta + beta,inf),'|| g*log(vr) - ci + cr - N''*lambda - C''*nabla - zeta + beta ||_inf');
    fprintf('%8.2g %s\n',norm(   C'*nabla + zeta + alpha,inf), '||  + C''*nabla + zeta + alpha ||_inf');
    fprintf('%8.2g %s\n',norm( - C'*nabla - zeta + beta,inf),'||  - C''*nabla - zeta + beta ||_inf');
else
    fprintf('%8.2g %s\n',norm(g.*reallog(vf) + ci + cf + N'*lambda + zeta - alpha,inf), '|| g*log(vf) + ci + cf + N''*lambda + zeta - alpha ||_inf');
    fprintf('%8.2g %s\n',norm(g.*reallog(vr) - ci + cr - N'*lambda - zeta - beta,inf),'|| g*log(vr) - ci + cr - N''*lambda - zeta - beta ||_inf');
end

fprintf('\n%s\n','Derived biochemical optimality conditions (fluxes)')

if isfield(model,'C')
    res = g.*reallog(vr./vf) + cr - cf - 2*ci - 2*N'*lambda - 2*C'*nabla - 2*zeta + beta - alpha;
    fprintf('%8.2g %s\n',norm(res,inf),'|| g*log(vr/vf) + cr - cf - 2*ci - 2*N''*lambda - 2*C''*nabla - 2*zeta + beta - alpha ||_inf');
else
    res = g.*reallog(vr./vf) + cr - cf - 2*ci - 2*N'*lambda - 2*zeta + beta + alpha;
    fprintf('%8.2g %s\n',norm(res,inf),'|| g*log(vr/vf) + cr - cf - 2*ci - 2*N''*lambda - 2*zeta + beta + alpha ||_inf');
end

if isfield(solution,'x')
    fprintf('\n%s\n','Optimality conditions (concentrations)')
    if detailed
        fprintf('%8.2g %s\n',norm(reallog(x),inf), '|| log(x)||_inf');
        fprintf('%8.2g %s\n',norm(u0,inf), '|| u0 ||_inf');
        fprintf('%8.2g %s\n',norm(lambda,inf), '|| lambda ||_inf');
        fprintf('%8.2g %s\n',norm(eta,inf), '|| eta ||_inf');
        fprintf('%8.2g %s\n',norm(delta,inf), '|| delta ||_inf');
        fprintf('%8.2g %s\n',norm(epsilon,inf), '|| epsilon ||_inf');
    end
    fprintf('%8.2g %s\n',norm(f.*reallog(x) + u0 - lambda + eta + delta,inf), '|| f.*log(x) + u0 - lambda + eta + delta ||_inf');
    fprintf('%8.2g %s\n',norm(f.*reallog(x0) + u0 + lambda - eta + epsilon,inf),'|| f.*log(x0) + u0 + lambda - eta + epsilon ||_inf');
    
    fprintf('\n%s\n','Derived biochemical optimality conditions (concentrations)')
    fprintf('%8.2g %s\n',norm(f.*reallog(x./x0) - 2*lambda + 2*eta + delta - epsilon,inf),'|| f.*log(x/x0) - 2*lambda + 2*eta + delta - epsilon ||_inf');
    fprintf('%8.2g %s\n',norm(f.*reallog(x.*x0) + 2*u0 + delta + epsilon,inf),'|| f.*log(x.*x0) + 2*u0 + delta + epsilon ||_inf');

    if isfield(model,'C')
        fprintf('\n%s\n','Derived biochemical optimality conditions (fluxes and concentrations)')
        if isfield(model,'C')
            fprintf('%8.2g %s\n',norm(g.*reallog(vf) + cf + ci + N'*(u0 + reallog(x) + eta + delta) + C'*nabla + zeta + alpha,inf),'|| g*log(vf) + cf + ci + N''*(u0 + f*log(x) + eta + delta) + C''*nabla + zeta + alpha ||_inf');
            fprintf('%8.2g %s\n',norm(g.*reallog(vr) + cr - ci - N'*(u0 + reallog(x) + eta + delta) - C'*nabla - zeta + beta,inf), '|| g*log(vr) + cr - ci - N''*(u0 + f*log(x) + eta + delta) - C''*nabla - zeta +  beta ||_inf');
        else
            fprintf('%8.2g %s\n',norm(g.*reallog(vf) + cf + ci + N'*(u0 + reallog(x) + eta + delta) + zeta + alpha,inf),'|| g*log(vf) + cf + ci + N''*(u0 + f*log(x) + eta + delta) + zeta + alpha ||_inf');
            fprintf('%8.2g %s\n',norm(g.*reallog(vr) + cr - ci - N'*(u0 + reallog(x) + eta + delta) - zeta + beta,inf),'|| g*log(vr) + cr - ci - N''*(u0 + f*log(x) + eta + delta) - zeta + beta ||_inf');
        end
    end
    
    if isfield(model,'C')
        fprintf('\n%s\n','Derived biochemical optimality conditions (fluxes and concentrations)')
        if isfield(model,'C')
            fprintf('%8.2g %s\n',norm(g.*reallog(vr./vf) + cr - cf - 2*ci - 2*N'*(u0 + f.*reallog(x) + eta + delta) - 2*C'*nabla - 2*zeta - alpha + beta,inf),   '|| g.*reallog(vr./vf) + cr - cf - 2*(ci + N''*(u0 + f*log(x)  + eta +   delta) + C''*nabla + zeta) - alpha + beta ||_inf');
            fprintf('%8.2g %s\n',norm(g.*reallog(vr./vf) + cr - cf - 2*ci + 2*N'*(u0 + f.*reallog(x0) - eta + epsilon) - 2*C'*nabla - 2*zeta - alpha + beta,inf),'|| g.*reallog(vr./vf) + cr - cf - 2*(ci - N''*(u0 + f*log(x0) - eta + epsilon) + C''*nabla + zeta) - alpha + beta ||_inf');
        else
            fprintf('%8.2g %s\n',norm(g.*reallog(vr./vf) + cr - cf - 2*ci - 2*N'*(u0 + f.*reallog(x) + eta + delta) - 2*zeta - alpha + beta,inf),   '|| g.*reallog(vr./vf) + cr - cf - 2*(ci + N''*(u0 + f*log(x)  + eta +   delta) + zeta) - alpha + beta ||_inf');
            fprintf('%8.2g %s\n',norm(g.*reallog(vr./vf) + cr - cf - 2*ci + 2*N'*(u0 + f.*reallog(x0) - eta + epsilon) - 2*zeta - alpha + beta,inf),'|| g.*reallog(vr./vf) + cr - cf - 2*(ci - N''*(u0 + f*log(x0) - eta + epsilon) + zeta) - alpha + beta ||_inf');
        end
    end
end

%%
fprintf('\n')
if isfield(solution,'v')
    vfpad = NaN*ones(size(model.S,2),1);
    vfpad(model.SConsistentRxnBool)=vf;
    vrpad = NaN*ones(size(model.S,2),1);
    vrpad(model.SConsistentRxnBool)=vr;
    T1 = table(model.lb, v,vfpad,vrpad,model.ub,'VariableNames',{'vl' 'v' 'vf' 'vr' 'vu'});
    if length(v)<=10
        disp(T1)
    end
end
%%
if isfield(solution,'x') 
    T2 = table(model.xl,x,model.xu,model.x0l,x0,model.x0u,delta,epsilon,...
        'VariableNames',{'xl' 'x' 'xu' 'x0l' 'x0' 'x0u' 'delta' 'epsilon'});
    T3 = table(model.dxl,dx,model.dxu,eta,...
        'VariableNames',{ 'dxl' 'dx' 'dxu' 'eta'});
    if length(x)<=10
        disp(T2)
        disp(T3)
    end
end

%%   
if 1
    %useful for debugging when tolerances not satifactory
    l = model.lb(model.SConsistentRxnBool);
    u = model.ub(model.SConsistentRxnBool);
    vi = vf - vr;
    logvrvf = reallog(vr./vf);

    figure
    histogram(log10(abs(zeta)))
    title('distribution of abs(zeta)')
    
    zBool = abs(zeta)>0.1;
    if any(l(zBool)~=0)
        figure
        histogram(l(zBool))
        title('distribution of internal lower bound when zeta substantial')
    end
    if any(u(zBool)~=inf)
        figure
        histogram(u(zBool))
        title('distribution of internal upper bound when zeta substantial')
    end
    
    if any(zBool)
        figure
        histogram(res(zBool))
        title('distribution of thermodynamic residual when zeta is substantial')
        
        figure
        plot(vi(zBool),logvrvf(zBool),'.')
        xlabel('Small net flux')
        ylabel('Change in chemical potential')
    end
    
    if any(~zBool)
        %thermo residual is large even when zeta is small
        figure
        histogram(res(~zBool))
        title('distribution of thermodynamic residual when zeta is small')
        
        
        figure
        plot(vi(~zBool),logvrvf(~zBool),'.')
        xlabel('Large net flux')
        ylabel('Change in chemical potential')
    end
    
    %%
    rBool = abs(res)>0.001;
    fwdBool = vi>0;
    if any(rBool)
        ab = alpha - beta;
        figure
        histogram(ab(rBool))
        title('distribution of alpha - beta when residual is large')
        
        %%
        %deviation from thermodynamic constraint occuring when net flux is
        %substantial
        figure
        res1 = g.*reallog(vr./vf) + cr - cf - 2*ci - 2*N'*lambda;
        histogram(res1(rBool))
        title('distribution of log(vr/vf) - 2*N''*lambda when residual is large')
        %%
        figure
        plot(vi(rBool),res(rBool),'.')
        xlabel('Net flux')
        ylabel('Large thermo residual')
        
        
        figure
        plot(vi(rBool & fwdBool),res(rBool & fwdBool),'.')
        xlabel('Net flux')
        ylabel('Large thermo residual')
        
        figure
        if isfield(model,'C')
            res5 = g.*reallog(vr./vf) + cr - cf - 2*ci - 2*N'*lambda - 2*C'*nabla - 2*zeta;
            histogram(res5(rBool))
            title('distribution of g.*log(vr/vf) + cr - cf - 2*ci - 2*N''*lambda - 2*C''*nabla - 2*zeta when residual is large')
        else
            res5 = g.*reallog(vr./vf) + cr - cf - 2*ci - 2*N'*lambda - 2*zeta;
            histogram(res5(rBool))
            title('distribution of g.*log(vr/vf) + cr - cf - 2*ci - 2*N''*lambda - 2*zeta when residual is large')
        end
    end
    
    if any(~rBool)
        figure
        plot(res(~rBool),vi(~rBool),'.')
        xlabel('Net flux')
        ylabel('Small thermo residual')
        
        figure
        plot(vi(~rBool & fwdBool),res(~rBool & fwdBool),'.')
        xlabel('Net flux')
        ylabel('Small thermo residual')
    end
    
    %%


end
