function taxon_id = i_species2taxid(species_name)

% ChatGPT code on 12/21/2022
%
%     % Replace spaces in the species name with underscores
%     species_name = strrep(species_name, ' ', '_');
% 
%     % Construct the API URL
%     api_url = sprintf('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=%s', species_name);
% 
%     % Retrieve the API response
%     api_response = webread(api_url);
% 
%     % Extract the Taxonomy ID from the API response
%     taxon_id = api_response.esearchresult.idlist{1};
% end


% function taxid = species2taxid(species)
    % Convert a species name to its NCBI Taxonomy ID
    %
    % Inputs:
    % - species: a string representing the species name
    %
    % Outputs:
    % - taxid: the NCBI Taxonomy ID for the species

    % Use webread to search the NCBI Taxonomy database for the species name
    url = sprintf('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=%s', species_name);
    result = webread(url);

    % Extract the NCBI Taxonomy ID from the search result
    taxid = regexp(result, '<Id>(\d+)<\/Id>', 'tokens');
    taxon_id = str2double(taxid{1});
end
