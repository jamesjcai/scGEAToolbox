function txt = i_stripcitation(txt)
%I_STRIPCITATION  Drop PMID citation tags from a label meant for display.
%
%   txt = gui.i_stripcitation(txt)
%
%   'UCell [PMID:34285779]'   ->  'UCell'
%   'Brennecke et al. (2013) [PMID:24056876]' -> 'Brennecke et al. (2013)'
%
% Method and menu names across the toolbox carry a [PMID:nnnn] tag so the
% source is discoverable where the choice is made. That is useful in a
% chooser and noise on an axis label, where the space is better spent on
% the quantity than on a citation the reader cannot click.
%
% Both bracket styles are handled, with or without a space after the colon.
% Nothing else is touched - descriptive parentheses such as
% 'AUCell (AUC recovery)' survive intact, since only PMID tokens match.
%
% Accepts char, string or cellstr and returns the same type.

if isempty(txt), return; end

wascell = iscell(txt);
waschar = ischar(txt);
s = string(txt);

s = regexprep(s, '\s*\[\s*PMID\s*:?\s*\d+\s*\]', '');
s = regexprep(s, '\s*\(\s*PMID\s*:?\s*\d+\s*\)', '');
s = strtrim(s);

if wascell
    txt = cellstr(s);
elseif waschar
    txt = char(s);
else
    txt = s;
end
end
