function main(inputData)

    % VARIABLES FOR INPUT PARAMS
    thresholdingApproach = 1;
    lowerGlobalThresholdVal = 0;
    upperGlobalThresholdVal = 0;
    trSheets = [];
    gmAndOperation = 'MIN';
    gmOrOperation = 'MAX';

    % RESULT HEADERS
    nonOptFilterLowlyHeaders1 = {'GeneId', 'Data', 'ExpressionLevel'};
    nonOptFilterLowlyHeaders2 = {'GeneId', 'Data', 'ExpressionLevel', 'ThApplied'};
    
    % TIMESTAMP
    t = datetime('now');
    strt = datestr(t);
    d = strrep(strt,':','-');
    
    % PREPARING INPUT PARAMS
    % (param) thresholdingApproach
    if inputData.localT1 == 1
        thresholdingApproach = 2;
    elseif inputData.localT2 == 1
        thresholdingApproach = 3;
    end
    % (param) lowerGlobalThresholdVal
    lowerGlobalThresholdVal = str2num(inputData.lowerGlobal);
    % (param) upperGlobalThresholdVal
    upperGlobalThresholdVal = str2num(inputData.upperGlobal);
    % (params) gene mapping AND/OR operators
    if inputData.gmGmSum == 1
        gmAndOperation = 'GM';
        gmOrOperation = 'SUM';
    elseif inputData.gmGmMax == 1
        gmAndOperation = 'GM';
    elseif inputData.gmMinSum == 1
        gmOrOperation = 'SUM';
    end
    % (params) boolean arrays for selected analysis tasks
    nonOptTasksBool = [inputData.filterHighlyLowlyExpressedGenes inputData.filterLowlyExpressedGenes inputData.ComparePhenotypeGenes];
    postOptTasksBool = [inputData.filterNonFluxReactions inputData.filterRateLimittingReactions inputData.calculateFluxShifts];
    
    
    % NON-OPT TASKS START (transcriptome analysis)
    if any(nonOptTasksBool(:) == 1)   
        
        % Get transcriptomics dataset sheet (phenotype) names
        trSheets = sheetnames(inputData.trDataPath);
        
        % Processing one sheet (phenotype) at a time
        if inputData.filterHighlyLowlyExpressedGenes == 1 || inputData.filterLowlyExpressedGenes == 1
            
            for i=1:1:height(trSheets) 
                % read a transcriptome dataset
                trData=readtable(inputData.trDataPath,'Sheet',trSheets{i});
                
                % percentile calculation (if selected)
                if inputData.percentile == 1        
                    lowerGlobalThresholdVal = calculatePercentile(trData.Data, str2double(inputData.lowerGlobal));
                    upperGlobalThresholdVal = calculatePercentile(trData.Data, str2double(inputData.upperGlobal));
                end
                
                % 1. GT1 (Global T1, no local rules)
                if inputData.globalT1 == 1 % 
         
                    % filterLowlyExpressedGenes GT1
                    if inputData.filterLowlyExpressedGenes == 1
                        lowlyExpressedGenes = findGenesBelowThresholdGT1(lowerGlobalThresholdVal, trData.Geneid, trData.Data);
                        folderName = strcat('resultsNonOptimization/lowlyExpressedGenes/',d,'/');
                        excelFileName = strcat(trSheets{i}, '_GT1_', string(lowerGlobalThresholdVal), '.xls');
                        saveResultTable(folderName, excelFileName, lowlyExpressedGenes, nonOptFilterLowlyHeaders1);
                    end
                    
                    % filterHighlyLowlyExpressedGenes GT1
                    if inputData.filterHighlyLowlyExpressedGenes == 1
                        geneExpressionLevels = findHighlyLowlyExpressedGenesGT1(lowerGlobalThresholdVal, trData.Geneid, trData.Data);
                        folderName = strcat('resultsNonOptimization/highlyLowlyExpressedGenes/',d,'/');
                        excelFileName = strcat(trSheets{i}, '_GT1_', string(lowerGlobalThresholdVal), '.xls');
                        saveResultTable(folderName, excelFileName, geneExpressionLevels, nonOptFilterLowlyHeaders1);
                    end
                    
                elseif inputData.localT1 == 1 % Must apply local rules (multiple tr datasets)
                    % filterLowlyExpressedGenes LT1
                    if inputData.filterLowlyExpressedGenes == 1
                        lowlyExpressedGenes = findGenesBelowThresholdLocal1(lowerGlobalThresholdVal, inputData.trDataPath, i);
                        folderName = strcat('resultsNonOptimization/lowlyExpressedGenes/',d,'/');
                        excelFileName = strcat(trSheets{i}, '_LT1_', string(lowerGlobalThresholdVal), '.xls');
                        saveResultTable(folderName, excelFileName, lowlyExpressedGenes, nonOptFilterLowlyHeaders2);
                    end
                    
                    % filterHighlyLowlyExpressedGenes LT1
                    if inputData.filterHighlyLowlyExpressedGenes == 1
                        geneExpressionLevels = findHighlyLowlyExpressedGenesLT1(lowerGlobalThresholdVal, inputData.trDataPath, i);
                        folderName = strcat('resultsNonOptimization/highlyLowlyExpressedGenes/',d,'/');
                        excelFileName = strcat(trSheets{i}, '_LT1_', string(lowerGlobalThresholdVal), '.xls');
                        saveResultTable(folderName, excelFileName, geneExpressionLevels, nonOptFilterLowlyHeaders2);
                    end

                else % inputData.localT2 == 1
                    % filterLowlyExpressedGenes LT2
                    if inputData.filterLowlyExpressedGenes == 1
                        lowlyExpressedGenes = findGenesBelowThresholdLocal2(lowerGlobalThresholdVal, upperGlobalThresholdVal, inputData.trDataPath, i);
                        folderName = strcat('resultsNonOptimization/lowlyExpressedGenes/',d,'/');
                        excelFileName = strcat(trSheets{i}, '_LT2_', string(lowerGlobalThresholdVal), '_', string(upperGlobalThresholdVal), '.xls');
                        saveResultTable(folderName, excelFileName, lowlyExpressedGenes, nonOptFilterLowlyHeaders2);
                    end
                    
                    % filterHighlyLowlyExpressedGenes LT2
                    if inputData.filterHighlyLowlyExpressedGenes == 1
                        geneExpressionLevels = findHighlyLowlyExpressedGenesLT2(lowerGlobalThresholdVal, upperGlobalThresholdVal, inputData.trDataPath, i);
                        folderName = strcat('resultsNonOptimization/highlyLowlyExpressedGenes/',d,'/');
                        excelFileName = strcat(trSheets{i}, '_LT2_', string(lowerGlobalThresholdVal), '_', string(upperGlobalThresholdVal), '.xls');
                        saveResultTable(folderName, excelFileName, geneExpressionLevels, nonOptFilterLowlyHeaders2);
                    end
                end
            end
        end        

        % Compare multiple sheets (phenotypes) at a time 
        if inputData.ComparePhenotypeGenes == 1
            % ComparePhenotypeGenes
            if string(inputData.nonOptCompareTarget) == "All"
                trSheets = sheetnames(inputData.trDataPath);
                for j=1:1:height(trSheets)
                    if trSheets{j} ~= string(inputData.nonOptCompareSource)
                        result = findUpDownRegulatedGenes(inputData.nonOptCompareSource, trSheets{j}, inputData.trDataPath);
                        folderName = strcat('resultsNonOptimization/geneExpressionLevelComparison/',d,'/');
                        excelFileName = strcat(inputData.nonOptCompareSource, '_compared_to_', trSheets{j}, '.xls');
                        saveResultTable(folderName, excelFileName, result, {'GeneId',strcat(inputData.nonOptCompareSource,'_Data(target)'),strcat(trSheets{j},'_Data(source)'),'Up/Down regulated'});
                    end
                end
            else
                result = findUpDownRegulatedGenes(inputData.nonOptCompareSource, inputData.nonOptCompareTarget, inputData.trDataPath);
                folderName = strcat('resultsNonOptimization/geneExpressionLevelComparison/',d,'/');
                excelFileName = strcat(inputData.nonOptCompareSource, '_compared_to_', inputData.nonOptCompareTarget, '.xls');
                saveResultTable(folderName, excelFileName, result, {'GeneId',strcat(inputData.nonOptCompareSource,'_Data(target)'),strcat(trSheets{j},'_Data(source)'),'Up/Down regulated'});
            end
        end
    end

    % Post-opt tasks
    if any(postOptTasksBool(:) == 1)
        
        % Start CobraToolbox
        if inputData.initCobraNoUpdates == 1
            initCobraToolbox(false);
        else
            initCobraToolbox(true);
        end
        
        % Integrate transcriptomics in the model
        if inputData.useExistingModels ~= 1
            createContextSpecificModel(inputData.modelPath,inputData.trDataPath,inputData.mediumDataPath,inputData.growthNotAffectingGeneDel, inputData.meetMinimumGrowthReq, thresholdingApproach, lowerGlobalThresholdVal, upperGlobalThresholdVal, inputData.objectiveFunction, gmAndOperation, gmOrOperation, inputData.constrAll, inputData.excludeBiomassEquation, inputData.biomassId, inputData.percentile);
        end
        
        if inputData.filterNonFluxReactions == 1
            dest = strcat('resultsPostOptimization\contextSpecificModels\*.xls');
            S = dir(dest);
            count = length(S);
            phenotypes = {};
            cnt = 1;
            
            for i=1:1:count
                temp = split(S(i).name,'_');
                phenotype = replace(temp(length(temp)),'.xls','');
                if ~any(strcmp(phenotypes,phenotype))
                    phenotypes(cnt) = phenotype;
                    filterNonFluxReactions(phenotype);
                    cnt = cnt + 1;
                end
            end
        end
        
        if inputData.filterRateLimittingReactions == 1
            dest = strcat('resultsPostOptimization\contextSpecificModels\*.xls');
            S = dir(dest);
            count = length(S);
            phenotypes = {};
            cnt = 1;
            
            for i=1:1:count
                temp = split(S(i).name,'_');
                phenotype = replace(temp(length(temp)),'.xls','');
                if ~any(strcmp(phenotypes,phenotype))
                    phenotypes(cnt) = phenotype;
                    cnt = cnt + 1;
                    filterRateLimittingReactions(phenotype);
                end
            end
        end
        
        if inputData.calculateFluxShifts == 1
            calculateFluxShifts(inputData.fluxShiftsSource, inputData.fluxShiftsTarget);
        end
       
    end
    f = msgbox('Operation Completed!');  
end