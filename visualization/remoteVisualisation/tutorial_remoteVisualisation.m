%% *Metabolic visualisation in ReconMap (Minerva) *
%% Authors: Alberto Noronha, Ines Thiele and Ronan M.T. Fleming
%% *Affiliation: *Luxembourg Centre for Systems Biomedicine, University of Luxembourg, Luxembourg.
%% Reviewer(s): Sylvain Arreckx, Alberto Nroonh
%% INTRODUCTION
% The visualisation of metabolic pathways is an essential tool to understand 
% the biological meaning underlying COBRA metabolic models. This would allow the 
% user to visualise what can not be appreciated at first sight by directly looking 
% to the model outputs. Here we present a visualisation through ReconMap [<#LyXCite-noronha_reconmap__2017 
% noronha_reconmap:_2017>], a virtual visualisation of human metabolism derived 
% from Recon 2.04 [<#LyXCite-thiele_community_driven_2013 thiele_community-driven_2013>]. 
% Diverse models and maps can be found in the website  <http://www.vmh.life http://www.vmh.life>.
%% *EQUIPMENT SETUP*
% In order to access remotely to ReconMap, the user has to be registered. To 
% obtain your credentials, you must access the ADMIN area and request an account.  
% To access ReconMap follow the link:  <http://www.vmh.life http://www.vmh.life> 
% .
% 
%      
% 
% Then, use your credentials to remotely access to <http://www.vmh.life http://www.vmh.life> 
%%
load('minerva.mat')
minerva.minervaURL = 'http://www.vmh.life/minerva/galaxy.xhtml';
minerva.map = 'ReconMap-2.01';
minerva.login = 'user_name';
minerva.password = 'user_password';
minerva.googleLicenseConsent = 'true';
%% 
% Initialise the Cobra Toolbox.
% 
% A specific solver might be required (depending on the analysis you want 
% to realise in the COBRA model).
%%
changeCobraSolver('gurobi','QP');
changeCobraSolver('gurobi', 'LP');
%% 
% Load your generic metabolic model. Recon's most recent version "Recon2.04" 
% can be freely downloaded from <http://www.vmh.life http://www.vmh.life>  

model = readCbModel('Recon2.v04.mat')
%% *PROCEDURE*
% *1. Overlay a flux distribution*
% 
% As an example of layout, we would like to see the the fluxes when maximizing 
% ATP production  through complex V (ATP synthase) in the Electron Transport Chain. 
% To do so, we use Flux Balance Analysis (FBA) and set as an objective function 
% the reaction responsible of this process (|'ATPS4m'). |
% 
% |ChangeObjective |function, changes the objective function of a constraint-based 
% model 
% 
% |optimizeCbModel |function solves a flux balance analysis problem.
%%
formula = printRxnFormula(model, 'ATPS4m')

model_atp_production = model  % re-name the model to do not modify the original one.
model_atp_production = changeObjective(model_atp_production, 'ATPS4m'); 
solution_atp_prod_max_regularised = optimizeCbModel(model_atp_production, 'max',1e-6);
solution_atp_prod_max_sparse = optimizeCbModel(model_atp_production, 'max','zero');
%% 
% This reaction is expressed in the regularised but not in the sparse, yet 
% both show up in the map!!!!!!!!!!!!!!!!!!!!!!

nnz(solution_atp_prod_max_regularised.v)
solution_atp_prod_max_regularised.v(strcmp(model.rxns,'r0603'))

nnz(solution_atp_prod_max_sparse.v)
solution_atp_prod_max_sparse.v(strcmp(model.rxns,'r0603'))
%% ANTICIPATED RESULTS
% |The buildFluxDistLayout |function, create a layout that is automatically 
% sent to the ReconMap website. After this, you can visualise your layout in <http://www.vmh.life 
% http://www.vmh.life>. Use your credentials to log in as it is previously explained.  
% Select your input map (minerva.model) and go to "overlays" section to find your 
% layout.   
% 
% 
%%
serverResponse = buildFluxDistLayout(minerva, model, solution_atp_production_max_regularised, 'atp_prod_max_regularised3',[],'#6617B5')
serverResponse = buildFluxDistLayout(minerva, model, solution_atp_production_max_sparse, 'atp_prod_max_sparse4',[],'#6617B5')
%% 
% If everything is correctly defined you should get a structure with 2 values. 
% If everything works fine, the output of this function should be:
% 
% |[1] 'Overlay was successfully sent to ReconMap!' |
% 
%  If there is any error, the message obtained will display:
% 
% |[0] 'Explanation of the obtained error' |
% 
% 
% 
% Note: If the "layout name" selected has been already given, an error might 
% appear. Please, try to give a new layout name each time you run the code. 
% 
% |[0]    'ERROR. Layout with given identifier ("atp_productio…'|
% 
% *2. Overlay a SubSystem*
% 
% There is also the possibility to highlight a specific subSystems by using 
% the function |generateSubsytemsLayout. A| subSystem is a group of metabolic 
% reactions involved in the same metabolic pathway, such as glycolysis, Oxidative 
% phosphorylation, citric acid cycle, etc.  Add the name of a specific subSystem 
% you want to highlight from the COBRA model (see the example, TCA cycle), and 
% the color reference. 
%%
generateSubsytemsLayout(minerva, model, 'Citric acid cycle', '#6617B5');
%% 
% 
% 
% Alternatively, the user can generate a layout of all common subSystems 
% between model and map using the function  |generateSubsystemLayouts. |
% 
% Note: every single layout can be observed individually, or merged with 
% other layouts. Therefore, making possible the visualization of several layouts 
% at the same time.
% 
% 
% 
% *RECONMAPS*
% 
% There is an additional set of maps available in VMH consistent with the 
% content of Recon 3D [3]. ReconMap3 is the general map that follows the same 
% approach as the previous iteration, but we have also included 6 organelle-specific 
% maps. You can see each of these by clicking on the corresponding button on the 
% interface as shown below.
% 
% 
% 
% To submit flux distributions to specific maps, users just need to change 
% the map variable from the minerva struct with the identifier of the desired 
% map using one of the lines of code displayed below:

minerva.map = 'ReconMap-3'; % ReconMap 3
minerva.map = 'reticulum'; % Endoplasmic reticulim map
minerva.map = 'peroxisome'; % Peroxisome map
minerva.map = 'nucleus'; % Nucleus map
minerva.map = 'mitochondrion'; % Mitochondrion map
minerva.map = 'lysosome'; % Lysosome map
minerva.map = 'golgi'; % Golgi apparatus map
%% REFERENCES: 
% [1] Alberto Noronha, Anna Dröfn Daníelsdóttir, Piotr Gawron, Freyr Jóhannsson, 
% Sofía Jónsdóttir, Sindri Jarlsson, Jón Pétur Gunnarsson, Sigurður Brynjólfsson, 
% Reinhard Schneider, Ines Thiele, and Ronan M. T. Fleming. ReconMap: an interactive 
% visualization of human metabolism. Bioinformatics , 33(4):605607, February 2017.
% 
% [2] Ines Thiele, Neil Swainston, Ronan M. T. Fleming, Andreas Hoppe, Swagatika 
% Sahoo, Maike K. Aurich, Hulda Haraldsdottir, Monica L. Mo, Ottar Rolfsson, Miranda 
% D. Stobbe, Stefan G. Thorleifsson, Rasmus Agren, Christian Bölling, Sergio Bordel, 
% Arvind K. Chavali, Paul Dobson, Warwick B. Dunn, Lukas Endler, David Hala, Michael 
% Hucka, Duncan Hull, Daniel Jameson, Neema Jamshidi, Jon J. Jonsson, Nick Juty, 
% Sarah Keating, Intawat Nookaew, Nicolas Le Novère, Naglis Malys, Alexander Mazein, 
% Jason A. Papin, Nathan D. Price, Evgeni Selkov Sr, Martin I. Sigurdsson, Evangelos 
% Simeonidis, Nikolaus Sonnenschein, Kieran Smallbone, Anatoly Sorokin, Johannes 
% H. G. M. van Beek, Dieter Weichart, Igor Goryanin, Jens Nielsen, Hans V. Westerho, 
% Douglas B. Kell, Pedro Mendes, and Bernhard Ø Palsson. A community-driven global 
% reconstruction of human metabolism. 5 Nat Biotech , 31(5):419425, May 2013.
% 
% 
% 
% [3] Elizabeth Brunk, Swagatika Sahoo, Daniel C Zielinski, Ali Altunkaya, 
% Andreas Dräger, Nathan Mih, Francesco Gatto, Avlant Nilsson, German Andres Preciat 
% Gonzalez, Maike Kathrin Aurich, Andreas Prlić, Anand Sastry, Anna D Danielsdottir, 
% Almut Heinken, Alberto Noronha, Peter W Rose, Stephen K Burley, Ronan M T Fleming, 
% Jens Nielsen, Ines Thiele & Bernhard O Palsson. Recon3D enables a three-dimensional 
% view of gene variation in human metabolism. _Nature Biotechnology_ volume 36, 
% pages 272–281 (2018)