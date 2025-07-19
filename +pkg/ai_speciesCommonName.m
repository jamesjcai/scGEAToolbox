function common = ai_speciesCommonName(scientific)
% speciesCommonName  Return the common name for a species
%
%   common = speciesCommonName(scientific)
%
%   Example:
%     speciesCommonName("Homo sapiens")
%     → "Human"
%
arguments
    scientific (1,:) char
end

% Normalize formatting
sci = lower(strtrim(scientific));

% Lookup table of known species
switch sci
    case "homo sapiens"
        common = "Human";
    case "mus musculus"
        common = "Mouse";
    case "danio rerio"
        common = "Zebrafish";
    case "drosophila melanogaster"
        common = "Fruit fly";
    % Add additional species mappings here

    otherwise
        % Return capitalized fallback of genus if unknown
        tokens = split(scientific);
        common = capitalize(tokens(1));
end

    function out = capitalize(txt)
        if isempty(txt)
            out = "";
        else
            out = [upper(txt(1)), lower(txt(2:end))];
        end
    end
end
