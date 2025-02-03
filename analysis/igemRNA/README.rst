IgemRNA
---------
*IgemRNA* is transcriptome analysis `MATLAB <https://se.mathworks.com/products/matlab.html?s_tid=hp_products_matlab>`_ software with a graphical user interface (GUI) designed   for the analysis of transcriptome data directly or the analysis of context-specific models generated from the provided model reconstruction, transciptome and optional medium     composition data files. *IgemRNA* facilitates some of the `Cobra Toolbox 3.0 <https://github.com/opencobra/cobratoolbox/>`_ constraint-based modelling functionalities for          context-specific model generation and performing optimisation methods like `FBA <https://opencobra.github.io/cobratoolbox/latest/modules/analysis/FBA/index.html>`_ 
or `FVA <https://opencobra.github.io/cobratoolbox/stable/modules/analysis/FVA/index.html>`_ on them.
Furthermore, *IgemRNA* can be used to validate transcriptomics data taking into account the interconnectivity 
of biochemical networks, steady state assumptions and Gene-Protein-Reaction (GPR) relationship. The result context-specific models can then be further analysed and compared to other phenotypes. The main advantage of *IgemRNA* software is that no previos programming skills are required due to the integrated user-friendly GUI which allows to select data files and data processing options in the *IgemRNA* form (see fig. 1). 

.. image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/IgemRNAForm.png
  :width: 550
**fig. 1** - Full IgemRNA form

Folder structure description
**********
The *IgemRNA* tool initially consists of 2 root folders (`Data <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/tree/main/Data>`_, `Scripts <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/tree/main/Scripts>`_), additional result folders are created when specific analysis tasks have been performed (*Results non-optimization, Results post-optimization*)   
and an `IgemRNA.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/IgemRNA.m>`_ file which calls the user graphical interface form. 
*Data* folder stores all input data files used for this tutorial and *Scripts* folder consists of all the script files that are being executed by the *IgemRNA* software according user-selected options in the *IgemRNA* form as well as the test cases provided in this demonstration.
The *Results non-optimization* (results of direct transcriptome analysis) and *Results post-optimization* (results of context-specific model analysis) are folders where all the result files are saved. 


Running IgemRNA software
**********
In order to run *IgemRNA* the user must first have the `MATLAB <https://se.mathworks.com/products/matlab.html?s_tid=hp_products_matlab>`_ environment installed and started as well as have the `IgemRNA <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA>`_ software downloaded and the files extracted. Then the user can navigate to the root folder of *IgemRNA* in *MATLAB* where the `IgemRNA.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/IgemRNA.m>`_ file is located and run it. 


Input data files
**********
To access all options in the *IgemRNA* form, the user must supply input data files. This can be done by pressing the *Open* button in the corresponding file upload row and selecting the files via *File Explorer*. **Transcriptomics data file** (see fig. 2B) is required to run *non-optimization tasks* but an additional **model reconstruction file** is required to access the *post-optimization tasks*. **Medium composition data file** (see fig. 2A) is optional, the selection of this data file does not extend the form, but specifies the provided exchange reaction constraints (upper and lower bounds) in the model for *post-optimization tasks* analysis.

Transcriptomics data and medium composition data can be provided in *.xls* or *.xlsx* formats, where columns are named respectively (see fig. 2) and sheet names correspond to dataset and phenotype names (*dataSetName_phenotypeName*). The model reconstruction file can be provided in *.xls*, *.xml* or other formats supported by *Cobra Toolbox 3.0.* **In case the transcriptomics data consists of multiple datasets, the genes listed in each dataset must be in the same order!**

.. image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/inputDataFormat.png
  :width: 500
**fig. 2** - Input data file structure; A - Medium data file structure; B - Transcriptomics data file structure

Data files used for this tutorial can be found in the `Data <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/tree/main/Data>`_ folder:

* `MediumData.xlsx <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA_v4/blob/main/Data/MediumData.xlsx>`_ (medium composition data)
* `Yeast_8_4_0.xls <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA_v4/blob/main/Data/Yeast_8_4_0.xls>`_ (the yeast consensus genome-scale model reconstruction)  
* `TranscriptomicsData.xlsx <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA_v4/blob/main/Data/TranscriptomicsData.xlsx>`_ (RNA-seq measurements), source `here <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130549>`_.

Non-optimization tasks 
**********

Non-optimization tasks include several transcriptomics data analysis tasks: 

* filter highly and lowly expressed genes, 
* filter lowly expressed genes, 
* filter up/down regulated genes between different phenotypes or data sets. 

The results of each task are stored in a different folder within the *Results non-optimization* folder: *Gene expression level comparison*, *Highly-lowly expressed genes*, *Lowly expressed genes* (see fig. 3).

.. image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/nonOptTasksFolderStructure.png
  :width: 400
**fig. 3** - *Non-optimization results* folder

**********
1. Filter highly and lowly expressed genes
**********
Non-optimization task *Filter highly and lowly expressed genes* generates result *Excel* files for each provided transcriptomics dataset. File names are assigned based on the provided dataset and phenotype names (from transcriptomics data), the applied thresholding approach (*GT1, LT1, LT2*, more on thresholding: `ThresholdingGeneMappingOptions.docx <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Supplementary%20Materials%202/ThresholdingGeneMappingOptions.docx>`_) and provided global thresholds values (see fig. 4). Each result *Excel* file contains one sheet with the list of genes provided by transcriptomics data and 4 columns: *GeneId*, *Data* (expression value), *ExpressionLevel* and *ThApplied*. The *ExpressionLevel* column contains the expression levels determined based on the selected thresholding approach, provided global and for thresholding approaches *LT1* and *LT2* calculated local thresholds. Column *ThApplied* displays whether a local or a global threshold for a specific gene was applied (see fig. 5). 

.. |pic1| image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/filterHighlyLowlyExpressedGenesFolder.PNG
   :width: 440

.. |pic2| image::  https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/filterHighlyLowlyExpressedGenesFile.PNG
   :width: 250

.. list-table::
   :widths: 200 200
   
   * - |pic1|
     - |pic2|

   * - **fig. 4** - *Filter highly and lowly expressed genes* folder
     - **fig. 5** - *Filter highly and lowly expressed genes* result file

There are two ways to perform this test case:

**1. Using GUI** - upload transciptomics data file, select a thresholding approach, input threshold value/s, select non-optimization tasks option *Filter highly and lowly expressed genes* and press *OK*.

**2. Run test case files** from the *Scripts* folder of *IgemRNA* tool:

* `TestCase_findHighlyLowlyExpressedGenesGT1.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Scripts/TestCase_findHighlyLowlyExpressedGenesGT1.m>`_;
* `TestCase_findHighlyLowlyExpressedGenesLT1.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Scripts/TestCase_findHighlyLowlyExpressedGenesLT1.m>`_;
* `TestCase_findHighlyLowlyExpressedGenesLT2.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Scripts/TestCase_findHighlyLowlyExpressedGenesLT2.m>`_.

**********
2. Filter lowly expressed genes 
**********
Non-optimization task *Filter lowly expressed genes* generates separate *Excel* result files for each dataset provided in transcriptomics data file. These result files contain filtered gene lists including genes with expression value below the given threshold value/s based on the selected thresholding approach. File names include dataset and phenotype names (from transcriptomics data file), the applied thresholding approach (*GT1, LT1, LT2*, more on thresholding: `ThresholdingGeneMappingOptions.docx <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Supplementary%20Materials%202/ThresholdingGeneMappingOptions.docx>`_) and provided global threshold value/s (see fig. 6). Each result file consists of 4 columns *GeneId*, *Data* (expression value), *ExpressionLevel* (in this case, only *Low*) and *ThApplied* to show whether a global or a local threshold was applied for a specific gene. Only those genes with expression values below the given threshold (depending on which thresholding approach is applied) are listed in the result files (see fig. 7).

.. |pic3| image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/filterLowlyExpressedGenesFolder.png
   :width: 400

.. |pic4| image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/filterLowlyExpressedGenesFile.PNG
   :width: 300
    
     
.. list-table::
   :widths: 200 200
   
   * - |pic3|
     - |pic4|

   * - **fig. 6** - *Filter lowly expressed genes* folder
     - **fig. 7** - *Filter lowly expressed genes* result file
     
     
There are two ways to perform this test case:

**1. Using GUI** - upload transciptomics data file, select a thresholding approach, input threshold value/s, select non-optimization tasks option *Filter lowly expressed genes* and press *OK*.

**2. Run test case files** from the *Scripts* folder of *IgemRNA* tool:

* `TestCase_findGenesBelowThresholdGT1.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Scripts/TestCase_findGenesBelowThresholdGT1.m>`_;
* `TestCase_findGenesBelowThresholdLocal1.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Scripts/TestCase_findGenesBelowThresholdLocal1.m>`_;
* `TestCase_findGenesBelowThresholdLocal2.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Scripts/TestCase_findGenesBelowThresholdLocal2.m>`_.

**********
3. Filter up/down regulated genes between phenotypes
**********
Non-optimization task *Filter up/down regulated genes between phenotypes* generates result *Excel* files in the *Gene expression level comparison* folder. Result file name contains dataset and phenotype names for both transcriptomics datasets that have been compared (see fig. 8). Result *Excel* data files contain a full gene list from the target dataset, expression values for both target and source datasets as well as the determined *up/down* regulation status (see fig. 9). 

.. |pic5| image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/filterUpDownRegulatedGenesBetweenPhenotypesFolder.png
   :width: 440

.. |pic6| image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/filterUpDownRegulatedGenesBetweenPhenotypesFile.png
   :width: 500
    
     
.. list-table::
   :widths: 200 200
   
   * - |pic5|
     - |pic6|

   * - **fig. 8** - *Filter up/down regulated genes between phenotypes* folder
     - **fig. 9** - *Filter up/down regulated genes between phenotypes* result file
     
There are two ways to perform this test case:

**1. Using GUI** - upload transciptomics data file with multiple datasets, select a thresholding approach, input threshold value/s, select non-optimization tasks option *Filter up/down regulated genes between phenotypes*, choose phenotypes to compare, and press *OK*.

**2. Run test case file** from the *Scripts* folder of *IgemRNA* tool: `TestCase_findUpDownRegulatedGenes.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Scripts/TestCase_findUpDownRegulatedGenes.m>`_.

Post-optimization tasks
**********
Context-specific models generated by *IgemRNA* post-optimization tasks as well as the results of the analysis performed on these models are saved in the *Results post-optimization* folder (see fig. 10). The results of post-optimization tasks are saved in folders with corresponding names: *Flux-shifts*, *Non-flux reactions* and *Rate limiting reactions*. 

.. image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/postOptTasksFolder.PNG
   :width: 400

**fig. 10** - *Post-optimization results* folder

There are two ways to generate models used in this tutorial:

**1. Using GUI** - upload transciptomics data file, model reconstrucion file and optionally a medium composition file, select a thresholding approach, input threshold value/s, select gene mapping approach and constraining options, choose one or more *post-optimization tasks* options and press *OK*.

**2. Run test case file** from the *Scripts* folder of *IgemRNA* tool: `TestCase_createContextSpecificModel.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Scripts/TestCase_generateContextSpecificModels.m>`_
 
  Since model generation and optimisation takes some time, especially for multiple transcriptomics datasets, context-specific model files used for this demonstration have  already been provided in the `Results post-optimization/Context-specific models <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/tree/main/Results%20post-optimization/Context-specific%20models>`_ folder.
  
**********
1. Filter non-flux reactions 
**********
Post-optimization task *Filter non-flux reactions* performs an analysis on created context-specific models of the same phenotype, the name of the phenotype is included in the result file name (see fig. 11). Each result *Excel* file contains a list of reactions that carry no flux in the result context-specific model (lower and upper bound equals 0). An additional sheet for all the common non-flux reactions of the same phenotype is also provided in the sheet *Common_(phenotypeName)* (see fig. 12).

.. |pic7| image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/nonFluxReactionsResultFolder.png
   :width: 450

.. |pic8| image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/nonFluxReactionTaskResultFile.png
   :width: 800
    
     
.. list-table::
   :widths: 455 805
   
   * - |pic7|
     - |pic8|

   * - **fig. 11** - *Filter non-flux reactions* folder
     - **fig. 12** - *Filter non-flux reactions* result file

There are two ways to perform this test case:

**1. Using GUI** - upload transciptomics data file, model reconstrucion file and optionally a medium composition file, select a thresholding approach, input threshold value/s, select gene mapping approach and constraining options and choose post-optimization task *Filter non-flux reactions* options. If context-specific models have already been generated, choose the *Use existing context-specific models* option. Press *OK*. More on thresholding and gene mapping: `ThresholdingGeneMappingOptions.docx <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Supplementary%20Materials%202/ThresholdingGeneMappingOptions.docx>`_.

**2. Run test case file** from the *Scripts* folder of *IgemRNA* tool: `TestCase_filterNonFluxReactions.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Scripts/TestCase_filterNonFluxReactions.m>`_
  
**********
2. Filter rate limiting reactions 
**********	
Post-optimization task *Filter rate limiting reactions* performs analysis on the generated context-specific models of the same phenotype, the phenotype name is included in the result files (see fig. 13). 

.. image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/filterRateLimitingReactionsResultFolder.png
   :width: 400

**fig. 13** - *Filter rate limiting reactions* folder

An analysis of these context-specific models have been performed in order to filter reactions that have reached the maximal flux value (*MaxFlux* calculated by `FVA <https://opencobra.github.io/cobratoolbox/stable/modules/analysis/FVA/index.html>`_) based on the upper bound constraint set according to transcriptomics data and GPR associations. An additional sheet *Common_(phenotypeName)* for common rate limiting reactions has also been added to the result file where rate limiting reactions that are present in all context-specific models of the same phenotype (see fig. 14).

.. image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/filterRateLimitingReactionsResultFile.png
   :width: 1000

**fig. 14** - *Filter rate limiting reactions* result file

There are two ways to perform this test case:

**1. Using GUI** - upload transciptomics data file, model reconstrucion file and optionally a medium composition file, select a thresholding approach, input threshold value/s, select gene mapping approach and constraining options and choose post-optimization task *Filter rate limiting reactions* options. If context-specific models have already been generated, choose the *Use existing context-specific models* option. Press *OK*. More on thresholding and gene mapping: `ThresholdingGeneMappingOptions.docx <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Supplementary%20Materials%202/ThresholdingGeneMappingOptions.docx>`_.

**2. Run test case file** from the *Scripts* folder of *IgemRNA* tool: `TestCase_filterRateLimittingReactions.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Scripts/TestCase_filterRateLimittingReactions.m>`_.

**********	
3. Flux shifts calculation between different phenotypesz
**********	
Post-optimization task *Calculate flux shifts between phenotypes* compares flux values (calculated by `FVA <https://opencobra.github.io/cobratoolbox/stable/modules/analysis/FVA/index.html>`_ on the context-specific models) between two different phenotypes. In this demonstration flux shifts analysis task was performed on the *S47D* phenotype datasets using *GT1* thresholding approach with the lower global threshold value of 0, *phenotype SRR8994358_WT* was compared to the wild type dataset *SRR8994357_WT* of the same thresholding approach and threshold values (see fig. 15). 

.. image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/calculateFluxShiftsResultFolder.png
   :width: 400

**fig. 15** - *Calculate flux shifts between phenotypes* folder

Each result file contains a full reaction list that corresponds to the *Reaction List* sheet in the provided model reconstruction file as well as additional columns for the calculation results: *MinFlux* and *MaxFlux* values (phenotype that is compared, fig. 16 L, M columns), *MinFlux/MaxFlux(dataset name)_(phenotype name)* phenotype that is used for comparison (fig. 16 N, O columns) and the *MinFlux/MaxFluxRatio* between these two phenotypes (fig. 16. P, Q columns). 

.. image:: https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/img/calculateFluxShiftsResultFile.png
   :width: 1000

**fig. 16** - *Calculate flux shifts between phenotypes* result file

There are two ways to perform this test case:

**1. Using GUI** - upload transciptomics data file, model reconstrucion file and optionally a medium composition file, select a thresholding approach, input threshold value/s, select gene mapping approach and constraining options and choose post-optimization task *Calculate flux shifts between phenotypes* options. If context-specific models have already been generated, choose the *Use existing context-specific models* option. Press *OK*. More on thresholding and gene mapping: `ThresholdingGeneMappingOptions.docx <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Supplementary%20Materials%202/ThresholdingGeneMappingOptions.docx>`_.

**2. Run test case file** from the *Scripts* folder of *IgemRNA* tool: `TestCase_calculateFluxShifts.m <https://github.com/BigDataInSilicoBiologyGroup/IgemRNA/blob/main/Scripts/TestCase_calculateFluxShifts.m>`_.


How to cite the IgemRNA software
**********	
[1] Kristina Grausa, Ivars Mozga, Karlis Pleiko, Agris Pentjuss, **Integrative Gene Expression and Metabolic Analysis tool IgemRNA**, https://doi.org/10.1101/2021.08.02.454732
