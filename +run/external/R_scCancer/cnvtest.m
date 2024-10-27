% https://raw.githubusercontent.com/wguo-research/scCancer/refs/heads/master/R/cnvFunction.R


function cnvList = prepareCNV(exprData, geneManifest, cellAnnotation, refData, species, genome)
    % Set default values for optional arguments
    if nargin < 4, refData = []; end
    if nargin < 5, species = 'human'; end
    if nargin < 6, genome = 'hg19'; end

    % Load gene chromosome information
    if strcmp(species, 'human')
        if strcmp(genome, 'hg38')
            geneChr = readtable('gene-chr-hg38.txt', 'ReadVariableNames', false);
        elseif strcmp(genome, 'hg19')
            geneChr = readtable('gene-chr-hg19.txt', 'ReadVariableNames', false);
        else
            error('Error in prepareCNV: genome "%s" is not allowed.', genome);
        end
    elseif strcmp(species, 'mouse')
        if strcmp(genome, 'mm10')
            geneChr = readtable('gene-chr-mm10.txt', 'ReadVariableNames', false);
        else
            error('Error in prepareCNV: genome "%s" is not allowed.', genome);
        end
    else
        error('Error in prepareCNV: species "%s" is not allowed.', species);
    end

    % If no reference data is provided, load default reference data
    if isempty(refData)
        if strcmp(species, 'human')
            refData = readmatrix('cnvRef_Data-HM.mat');
        elseif strcmp(species, 'mouse')
            refData = readmatrix('cnvRef_Data-boneMarrow-MS.mat');
        end
    end

    % Create reference annotation
    refAnno = table(string(cellstr(colnames(refData))), repmat("Reference", size(refData, 2), 1), 'VariableNames', {'cellName', 'cellAnno'});

    % Combine expression data and reference data
    commonGenes = intersect(rownames(exprData), rownames(refData));
    exprData = exprData(commonGenes, :);
    refData = refData(commonGenes, :);
    combinedData = [exprData refData];

    % Update row names with Ensembl IDs
    combinedData.Properties.RowNames = geneManifest{commonGenes, 'EnsemblID'};

    % Create cell annotation with combined reference annotation
    cellAnno = [cellAnnotation; refAnno];
    cellAnno.Properties.RowNames = cellAnno.cellName;

    % Filter gene chromosome data based on common genes
    commonGenesChr = intersect(geneChr.EnsemblID, combinedData.Properties.RowNames);
    geneChr = geneChr(ismember(geneChr.EnsemblID, commonGenesChr), :);

    % Order by chromosome start position
    combinedData = combinedData(geneChr.EnsemblID, :);

    % Return structured data
    cnvList.exprData = combinedData;
    cnvList.geneChr = geneChr;
    cnvList.cellAnno = cellAnno;
end


function cnvList = rmGeneForCNV(cnvList, cutoff, minCell)
    % Calculate the mean and sum for gene selection
    geneMean = mean(cnvList.exprData, 2);
    geneSum = sum(cnvList.exprData > 0, 2);

    % Select genes based on criteria
    genesSel = (geneMean >= cutoff) & (geneSum >= minCell);

    % Filter expression data and chromosome info
    cnvList.exprData = cnvList.exprData(genesSel, :);
    cnvList.geneChr = cnvList.geneChr(genesSel, :);
end


