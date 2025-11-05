%% MicroMap exploration with CellDesigner
%% Author: Cyrille C. Thinnes, University of Galway
% The aim of this tutorial is to guide you on *how to inspect the MicroMap with 
% CellDesigner*. CellDesigner will require a fair bit of RAM to handle the MicroMap. 
% We therefore suggest you use a computer with at least 16Gb of RAM. You can find 
% an <https://youtu.be/K6YvZC_gD6c accompanying video walkthrough on YouTube>. 
%% Downloading resources
% You can download the MicroMap CellDesigner .xml file from the <https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FZKMJ8 
% Harvard Dataverse>. We also provide an image-based .pdf version, i.e., you can 
% explore it visually, but you cannot search species within. You may download 
% the latest CellDesigner version from the <https://celldesigner.org/ CellDesigner 
% website>, which also contains plenty of <https://celldesigner.org/documents.html 
% technical documentation>. In this tutorial, we shall focus on specifically exploring 
% the MicroMap.
%% Inspecting the MicroMap .pdf version
% For inspecting the .pdf, we suggest you right-click the file in your file 
% explorer and choose ‘open with’ your internet browser. We found this to be the 
% fastest way to explore the .pdf, as using a dedicated .pdf viewer, such as Adobe 
% Acrobat, was noticeably slower.
%% Inspecting the MicroMap .xml with CellDesigner
% General navigation
% Open CellDesigner, then open the MicroMap .xml from within CellDesigner, by 
% clicking on the ‘Open’ icon, or by using File>Open (*Figure 1*). Opening the 
% MicroMap directly from your file explorer, e.g., by double-clicking, does not 
% work – ensure you open the file from within CellDesigner. The MicroMap may take 
% a moment to load, as indicated by the blue animation bar in the right bottom 
% corner. The MicroMap is fully loaded when the animation stops and you can see 
% a list of species in the list view. Most likely, the first view will be zoomed 
% all the way into a white section at the top left corner.
% 
% To start exploring, you can grab one of the scroll handles and pan around 
% the map (*Figure 1*). We suggest grabbing the handles directly for moving around 
% the map – using your scroll wheel will likely be slower. You may zoom in and 
% out of the MicroMap by clicking on the dedicated icons at the top (*Figure 1*). 
% You will see an icon’s function by hovering over it.
% 
% 
% 
% *Figure 1:* Overview of CellDesigner control functions to inspect the MicroMap. 
% Open your .xml file from within CellDesigner. Use the Zoom icons to adjust your 
% zoom levels. Pan by grabbing the scroll handles.
% Metabolites
% You can search for metabolites by using the ‘find species’ function (control 
% + F, or Component > Find Species in the menu bar). A window will appear and 
% allow you to select different search options (*Figure 2*). Please note that 
% you will need to select ‘name’ to search for a VMH ID. CellDesigner IDs are 
% arbitrarily assigned by CellDesigner, i.e., a VMH ID will correspond to a CellDesigner 
% name, not CellDesigner ID! For example, for searching an exact match for cytosolic 
% oxoglutaric acid, you shall enter akg[c] into a species name search, and select 
% ‘equal’ in ‘pattern’. Your found metabolite will be highlighted with a pink 
% bounding box, with any associated reactions highlighted in blue. Please note 
% that several reactions can be linked to a single found metabolite. As a metabolite 
% may occur several times throughout the map, you can navigate different search 
% hits by clicking the ‘Next’ and ‘Prev’ buttons. The metabolite will be highlighted 
% at the current zoom level. Therefore, please ensure to be at a reasonably zoomed-in 
% level where you can easily spot the pink bounding box. If you still cannot spot 
% your found metabolite, it may be barely off screen – you may need to scroll 
% a wee bit to fully see it. You will get a ‘not matched’ popup window if your 
% desired metabolite was not found.
% 
% 
% 
% *Figure 2:* Overview of metabolite search. Open the ‘Find Species’ dialogue 
% through Ctrl + F and search for VMH IDs by ensuring you have ‘name’ ticked in 
% the search box. The found metabolites will be highlighted in a pink bounding 
% box, with their associated reactions highlighted in blue. Cycle through the 
% search hits by using the ‘Prev’ and ‘Next’ buttons.
% Reactions
% The most straightforward way to identify a given reaction within CellDesigner 
% is to click on its reaction arrow. The selected reaction information, including 
% its VMH ID (= CellDesigner name), will be displayed in the notes window, and 
% the highlighted reaction will point to all the associated metabolites. This 
% method is particularly helpful in areas where, due to CellDesigner layout constraints, 
% several reaction arrows are close to each other, overlapping, or which span 
% distantly distributed metabolites.
% 
% Reaction names are displayed next to their corresponding reaction arrow, where 
% possible. However, due to naming constraints within CellDesigner, this was not 
% possible for reaction IDs that start with a number, which, instead, show an 
% automatically generated reaction label. This limitation is cosmetic only, as 
% each reaction remains associated with its VMH ID, which will be displayed when 
% clicking on the reaction arrow, as outlined above.
% 
% Unlike metabolites, however, reactions are not searchable. For finding a specific 
% reaction, click on the ‘Reactions’ tab within the list window to display an 
% alphabetically ordered list of all contained VMH reaction IDs, which will then 
% allow you to find and select your reaction of interest (*Figure 3*). Clicking 
% on the reaction ID will jump to the relevant section in the MicroMap and highlight 
% the reaction of interest. As with searching for metabolites, the found reaction 
% will be highlighted at the current zoom level. Therefore, please ensure to be 
% at a reasonably zoomed-in display before clicking on your reaction of interest.
%% 
% *Figure 3:* Overview of reaction search. Click on the ‘Reactions’ tab in the 
% list view to access the alphabetically ordered reaction list. Select the reaction 
% of interest, which will then appear highlighted on the map.
% 
% We encourage you to play around with the MicroMap and CellDesigner functionality 
% to explore and identify the most suitable workflows for your use-case. For more 
% information on CellDesigner and potential troubleshooting, please consult the 
% ample <https://celldesigner.org/documents.html CellDesigner documentation> online.