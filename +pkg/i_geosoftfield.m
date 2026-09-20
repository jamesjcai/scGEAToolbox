function v = i_geosoftfield(lines, key)
% I_GEOSOFTFIELD  Read an attribute out of a GEO SOFT record.
%
%   v = pkg.i_geosoftfield(lines, key)
%
%   Inputs:
%     lines - SOFT record lines, as returned by pkg.i_geosoftrecord
%     key   - attribute name fragment, matched case-insensitively against
%             the attribute name with its entity prefix stripped. "organism"
%             therefore finds '!Sample_organism_ch1' in a sample record and
%             '!Series_sample_organism' in a series record.
%
%   Outputs:
%     v     - string scalar. Attributes that repeat (characteristics,
%             protocols, contributors) are de-duplicated and joined with
%             '; '. V is "" when no attribute matches.
%
%   See also pkg.i_geosoftrecord

    lines = string(lines);
    lines = lines(startsWith(lines, "!") & contains(lines, "="));
    if isempty(lines)
        v = "";
        return;
    end

    names = regexprep(strtrim(extractBefore(lines, "=")), '^!\w+?_', '');
    hit = contains(names, key, 'IgnoreCase', true);
    if ~any(hit)
        v = "";
        return;
    end

    v = strtrim(extractAfter(lines(hit), "="));
    v = strjoin(unique(v, 'stable'), "; ");
end
