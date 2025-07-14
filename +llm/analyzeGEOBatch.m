function results = analyzeGEOBatch(accessions, varargin)
    % ANALYZEGEOBATCH Batch analysis of multiple GEO records
    %
    % SYNTAX:
    %   results = analyzeGEOBatch(accessions)
    %   results = analyzeGEOBatch(accessions, 'SaveResults', true)
    %
    % INPUT:
    %   accessions - Cell array of GEO accession numbers
    %   
    % Optional Parameters:
    %   'SaveResults' - Save results to MAT file (default: false)
    %   'OutputFile' - Output filename (default: 'geo_analysis_results.mat')
    %   'Verbose' - Display progress (default: true)
    %
    % OUTPUT:
    %   results - Structure array containing analysis results for each accession
    %
    % EXAMPLE:
    %   accessions = {'GSM7855468', 'GSM7855469', 'GSM7855470'};
    %   results = llm.analyzeGEOBatch(accessions, 'SaveResults', true);
    
    % Parse input arguments
    p = inputParser;
    addRequired(p, 'accessions', @(x) iscell(x) || isstring(x));
    addParameter(p, 'SaveResults', false, @islogical);
    addParameter(p, 'OutputFile', 'geo_analysis_results.mat', @ischar);
    addParameter(p, 'Verbose', true, @islogical);
    parse(p, accessions, varargin{:});
    
    % Convert to cell array if needed
    if isstring(accessions)
        accessions = cellstr(accessions);
    end
    
    % Initialize results
    results = struct();
    
    % Process each accession
    for i = 1:length(accessions)
        if p.Results.Verbose
            fprintf('Processing %d/%d: %s\n', i, length(accessions), accessions{i});
        end
        
        try
            % Extract information for this accession
            geoInfo = extractGEOInfo(accessions{i});
            results(i).accession = accessions{i};
            results(i).data = geoInfo;
            results(i).success = true;
            results(i).timestamp = datetime('now');
            
        catch ME
            if p.Results.Verbose
                warning(ME.identifier, 'Failed to process %s: %s', accessions{i}, ME.message);
            end
            results(i).accession = accessions{i};
            results(i).data = struct();
            results(i).success = false;
            results(i).error = ME.message;
            results(i).timestamp = datetime('now');
        end
    end
    
    % Save results if requested
    if p.Results.SaveResults
        save(p.Results.OutputFile, 'results', 'accessions');
        if p.Results.Verbose
            fprintf('Results saved to %s\n', p.Results.OutputFile);
        end
    end
    
    % Generate batch summary
    if p.Results.Verbose
        generateBatchSummary(results);
    end
end

function generateBatchSummary(results)
    % Generate summary statistics for batch analysis
    fprintf('\n--- Batch Analysis Summary ---\n');
    fprintf('Total records processed: %d\n', length(results));
    fprintf('Successful analyses: %d\n', sum([results.success]));
    fprintf('Failed analyses: %d\n', sum(~[results.success]));
    
    % Extract research domains from successful analyses
    domains = {};
    organisms = {};
    experiment_types = {};
    
    for i = 1:length(results)
        if results(i).success && isfield(results(i).data, 'llm_analysis')
            if isfield(results(i).data.llm_analysis, 'research_domain')
                domains{end+1} = results(i).data.llm_analysis.research_domain;
            end
        end
        
        if results(i).success && isfield(results(i).data, 'structured')
            if isfield(results(i).data.structured, 'organism')
                organisms{end+1} = results(i).data.structured.organism;
            end
            if isfield(results(i).data.structured, 'library_strategy')
                experiment_types{end+1} = results(i).data.structured.library_strategy;
            end
        end
    end
    
    % Display unique domains and organisms
    if ~isempty(domains)
        fprintf('\nResearch domains identified:\n');
        uniqueDomains = unique(domains);
        for i = 1:length(uniqueDomains)
            fprintf('  - %s\n', uniqueDomains{i});
        end
    end
    
    if ~isempty(organisms)
        fprintf('\nOrganisms studied:\n');
        uniqueOrganisms = unique(organisms);
        for i = 1:length(uniqueOrganisms)
            fprintf('  - %s\n', uniqueOrganisms{i});
        end
    end
    
    if ~isempty(experiment_types)
        fprintf('\nExperiment types:\n');
        uniqueTypes = unique(experiment_types);
        for i = 1:length(uniqueTypes)
            fprintf('  - %s\n', uniqueTypes{i});
        end
    end
end

function displayGEOInfo(geoInfo)
    % DISPLAYGEOINFO Display extracted GEO information in a formatted way
    %
    % SYNTAX:
    %   displayGEOInfo(geoInfo)
    %
    % INPUT:
    %   geoInfo - Structure returned by extractGEOInfo
    
    if isempty(geoInfo) || ~isstruct(geoInfo)
        fprintf('No valid GEO information to display.\n');
        return;
    end
    
    fprintf('\n=== GEO Record Analysis ===\n');
    fprintf('Accession: %s\n', geoInfo.accession);
    fprintf('URL: %s\n', geoInfo.url);
    
    % Display structured information
    if isfield(geoInfo, 'structured')
        fprintf('\n--- Structured Information ---\n');
        s = geoInfo.structured;
        fields = fieldnames(s);
        for i = 1:length(fields)
            field = fields{i};
            value = s.(field);
            fprintf('%s: %s\n', upper(strrep(field, '_', ' ')), value);
        end
    end
    
    % Display LLM analysis
    if isfield(geoInfo, 'llm_analysis')
        fprintf('\n--- LLM Analysis ---\n');
        if isfield(geoInfo.llm_analysis, 'research_domain')
            fprintf('Research Domain: %s\n', geoInfo.llm_analysis.research_domain);
        end
        if isfield(geoInfo.llm_analysis, 'research_objectives')
            fprintf('\nResearch Objectives:\n%s\n', geoInfo.llm_analysis.research_objectives);
        end
        if isfield(geoInfo.llm_analysis, 'detailed_analysis')
            fprintf('\nDetailed Analysis:\n%s\n', geoInfo.llm_analysis.detailed_analysis);
        end
    end
    
    % Display text analysis
    if isfield(geoInfo, 'text_analysis')
        fprintf('\n--- Text Analysis ---\n');
        if isfield(geoInfo.text_analysis, 'top_keywords')
            fprintf('Top Keywords: %s\n', strjoin(geoInfo.text_analysis.top_keywords(1:min(10, end)), ', '));
        end
        if isfield(geoInfo.text_analysis, 'word_count')
            fprintf('Word Count: %d\n', geoInfo.text_analysis.word_count);
        end
    end
    
    % Display summary
    if isfield(geoInfo, 'summary')
        fprintf('\n--- Summary ---\n');
        s = geoInfo.summary;
        fields = fieldnames(s);
        for i = 1:length(fields)
            field = fields{i};
            value = s.(field);
            if iscell(value)
                fprintf('%s: %s\n', upper(strrep(field, '_', ' ')), strjoin(value, ', '));
            else
                fprintf('%s: %s\n', upper(strrep(field, '_', ' ')), string(value));
            end
        end
    end
end

function exportToTable(results, filename)
    % EXPORTTOTABLE Export batch analysis results to a table
    %
    % SYNTAX:
    %   exportToTable(results, filename)
    %
    % INPUT:
    %   results - Structure array from analyzeGEOBatch
    %   filename - Output filename (CSV or Excel)
    
    % Extract key information into table format
    accessions = {results.accession}';
    success = [results.success]';
    timestamps = [results.timestamp]';
    
    % Initialize other columns
    titles = cell(length(results), 1);
    organisms = cell(length(results), 1);
    experiment_types = cell(length(results), 1);
    research_domains = cell(length(results), 1);
    
    % Fill data from successful analyses
    for i = 1:length(results)
        if results(i).success && isfield(results(i).data, 'structured')
            s = results(i).data.structured;
            titles{i} = s.title;
            organisms{i} = s.organism;
            experiment_types{i} = s.library_strategy;
        end
        
        if results(i).success && isfield(results(i).data, 'llm_analysis')
            if isfield(results(i).data.llm_analysis, 'research_domain')
                research_domains{i} = results(i).data.llm_analysis.research_domain;
            end
        end
    end
    
    % Create table
    T = table(accessions, success, timestamps, titles, organisms, ...
              experiment_types, research_domains, ...
              'VariableNames', {'Accession', 'Success', 'Timestamp', 'Title', ...
                               'Organism', 'ExperimentType', 'ResearchDomain'});
    
    % Write to file
    [~, ~, ext] = fileparts(filename);
    if strcmpi(ext, '.csv')
        writetable(T, filename);
    else
        writetable(T, filename, 'FileType', 'spreadsheet');
    end
    
    fprintf('Results exported to %s\n', filename);
end

% Example usage function
function runExample()
    % RUNEXAMPLE Demonstrate the GEO analysis functions
    
    fprintf('=== GEO Analysis Example ===\n\n');
    
    % Single record analysis
    fprintf('1. Analyzing single GEO record...\n');
    geoInfo = llm.extractGEOInfo('GSM7855468');
    displayGEOInfo(geoInfo);
    
    % Batch analysis (uncomment to run)
    % fprintf('\n2. Batch analysis example...\n');
    % accessions = {'GSM7855468', 'GSM7855469'};
    % results = analyzeGEOBatch(accessions, 'Verbose', true);
    % exportToTable(results, 'geo_analysis_results.csv');
end