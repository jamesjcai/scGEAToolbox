function [T] = fetch_scRNAseq_GSEs

    % Change this to your real email
    email = 'jamescail@genomezoo.net';    
    % Search query
    query = '"single cell RNA-seq"[All Fields] OR scRNA-seq[All Fields] OR snRNA-seq[All Fields]';
    
    % Base URLs
    base_esearch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi';
    base_esummary = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi';
    
    % Step 1: ESearch to get list of GDS IDs (up to 10,000)
    retmax = 10;  % Increase if needed
    search_url = sprintf('%s?db=gds&term=%s&retmax=%d&email=%s', ...
        base_esearch, urlencode(query), retmax, email);
    
    disp('Searching GEO...');
    xml = webread(search_url);
    idlist = regexp(xml, '<Id>(\d+)</Id>', 'tokens');
    idlist = [idlist{:}];
    
    fprintf('Found %d records. Fetching details...\n', numel(idlist));
    
    % Step 2: For each ID, call ESummary to get accession and title
    results = strings(numel(idlist), 2);
    
    for i = 1:numel(idlist)
        id = idlist{i};
        summary_url = sprintf('%s?db=gds&id=%s&email=%s', base_esummary, id, email);
        pause(3);
        xml = webread(summary_url);
        
        acc = regexp(xml, '<Item Name="Accession" Type="String">([^<]+)</Item>', 'tokens', 'once');
        title = regexp(xml, '<Item Name="title" Type="String">([^<]+)</Item>', 'tokens', 'once');
        
        if ~isempty(acc)
            results(i,1) = acc{1};
        else
            results(i,1) = "";
        end
        if ~isempty(title)
            results(i,2) = title{1};
        else
            results(i,2) = "";
        end
    end
    
    % Step 3: Write to CSV
    T = table(results(:,1), results(:,2), 'VariableNames', {'Accession','Title'});
    %writetable(T, 'scRNAseq_GSEs.csv');    
    %fprintf('Saved %d GSE accession numbers to scRNAseq_GSEs.csv\n', size(T,1));
end

%{
function fetch_scRNAseq_GSEs(email)
%FETCH_SCRNASEQ_GSES Fetch all GSE accession numbers related to scRNA-seq from GEO (batched, faster)
%   Usage:
%       fetch_scRNAseq_GSEs('your.email@example.com')
%
%   The function searches GEO for common scRNA-seq keywords,
%   retrieves accession numbers and titles, and saves them
%   to 'scRNAseq_GSEs.csv' in the current folder.

% Define search query
query = '"single cell RNA-seq"[All Fields] OR scRNA-seq[All Fields] OR snRNA-seq[All Fields]';

% Base URLs
base_esearch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi';
base_esummary = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi';

% Step 1: ESearch to get list of GDS IDs (up to 10,000)
retmax = 10000;  % Increase if needed
search_url = sprintf('%s?db=gds&term=%s&retmax=%d&email=%s', ...
    base_esearch, urlencode(query), retmax, email);

disp('Searching GEO...');
xml = webread(search_url);
idlist = regexp(xml, '<Id>(\d+)</Id>', 'tokens');
idlist = [idlist{:}];

n = numel(idlist);
fprintf('Found %d records. Fetching details in batches...\n', n);

% Step 2: Fetch details in batches
batch_size = 100;  % can go up to ~200 safely
results = strings(n, 2);

for startIdx = 1:batch_size:n
    endIdx = min(startIdx+batch_size-1, n);
    batch_ids = strjoin(idlist(startIdx:endIdx), ',');
    
    summary_url = sprintf('%s?db=gds&id=%s&email=%s', base_esummary, batch_ids, email);
    xml = webread(summary_url);
    
    % Parse each record in batch
    accs = regexp(xml, '<Item Name="Accession" Type="String">([^<]+)</Item>', 'tokens');
    titles = regexp(xml, '<Item Name="title" Type="String">([^<]+)</Item>', 'tokens');
    
    % Sometimes fewer records returned; match carefully
    for j = 1:length(accs)
        idx = startIdx + j - 1;
        if idx <= n
            if ~isempty(accs{j})
                results(idx,1) = accs{j}{1};
            end
            if ~isempty(titles{j})
                results(idx,2) = titles{j}{1};
            end
        end
    end
end

% Step 3: Write to CSV
T = table(results(:,1), results(:,2), 'VariableNames', {'Accession','Title'});
writetable(T, 'scRNAseq_GSEs.csv');

fprintf('Done! Saved %d GSE accession numbers to scRNAseq_GSEs.csv\n', size(T,1));
end
%}


%{
function [T] = fetch_scRNAseq_GSEs
%FETCH_SCRNASEQ_GSES Fetch GSE accession numbers related to scRNA-seq from GEO (batched)
%   Saves accession, title, summary, and publication date to CSV.
%   Usage:
%       fetch_scRNAseq_GSEs('your.email@example.com')

email = 'jamescail@genomezoo.net'; 
% Define search query
query = '"single cell RNA-seq"[All Fields] OR scRNA-seq[All Fields] OR snRNA-seq[All Fields]';

% Base URLs
base_esearch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi';
base_esummary = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi';

% Step 1: ESearch to get list of GDS IDs (up to 10,000)
retmax = 10000;
search_url = sprintf('%s?db=gds&term=%s&retmax=%d&email=%s', ...
    base_esearch, urlencode(query), retmax, email);

disp('Searching GEO...');
xml = webread(search_url);
idlist = regexp(xml, '<Id>(\d+)</Id>', 'tokens');
idlist = [idlist{:}];

n = numel(idlist);
fprintf('Found %d records. Fetching details in batches...\n', n);

% Step 2: Fetch details in batches
batch_size = 10;
accessions = strings(n,1);
titles     = strings(n,1);
summaries  = strings(n,1);
pubdates   = strings(n,1);

for startIdx = 1:batch_size:n
    endIdx = min(startIdx+batch_size-1, n);
    batch_ids = strjoin(idlist(startIdx:endIdx), ',');
    
    summary_url = sprintf('%s?db=gds&id=%s&email=%s', base_esummary, batch_ids, email);
    xml = webread(summary_url);
    
    %assignin('base','xml',xml)
    % Parse each field
    accs = regexp(xml, '<Item Name="Accession" Type="String">([^<]+)</Item>', 'tokens');
    ts   = regexp(xml, '<Item Name="title" Type="String">([^<]+)</Item>', 'tokens');
    sums = regexp(xml, '<Item Name="summary" Type="String">([^<]+)</Item>', 'tokens');
    pubs = regexp(xml, '<Item Name="PDAT" Type="String">([^<]+)</Item>', 'tokens');


    % pause
    
    for j = 1:length(accs)
        idx = startIdx + j - 1;
        if idx <= n
            if ~isempty(accs{j})
                accessions(idx) = accs{j}{1};
            end
            if ~isempty(ts{j})
                titles(idx) = ts{j}{1};
            end
            if ~isempty(sums{j})
                summaries(idx) = sums{j}{1};
            end
            if ~isempty(pubs{j})
                pubdates(idx) = pubs{j}{1};
            end
        end
    end
end

% Step 3: Write to CSV
T = table(accessions, titles, summaries, pubdates, ...
    'VariableNames', {'Accession','Title','Summary','PubDate'});
%writetable(T, 'scRNAseq_GSEs.csv');

%fprintf('Done! Saved %d GSEs to scRNAseq_GSEs.csv\n', size(T,1));
end


%}


%{
from Bio import Entrez
import csv

# Always tell NCBI who you are
Entrez.email = "jamescail@genomezoo.net"  # <-- change this!

# Define search query for scRNA-seq studies
query = '"single cell RNA-seq"[All Fields] OR scRNA-seq[All Fields] OR snRNA-seq[All Fields] OR scRNAseq[All Fields]'

# Search in GEO DataSets (db="gds") for series (GSE)
print("Searching GEO...")
handle = Entrez.esearch(db="gds", term=query, retmax=10000)  # adjust retmax if you expect more
record = Entrez.read(handle)
id_list = record["IdList"]

print(f"Found {len(id_list)} records. Fetching details...")

results = []
for gds_id in id_list:
    summary_handle = Entrez.esummary(db="gds", id=gds_id)
    summary = Entrez.read(summary_handle)
    acc = summary[0]["Accession"]
    title = summary[0]["title"]
    taxon = summary[0]["taxon"]
    results.append((acc, title, taxon))

# Write to CSV
with open("scRNAseq_GSEs.csv", mode="w", newline="", encoding="utf-8") as file:
    writer = csv.writer(file)
    writer.writerow(["Accession", "Title"])
    writer.writerows(results)

print(f"Saved {len(results)} GSE accession numbers to scRNAseq_GSEs.csv")


%}


