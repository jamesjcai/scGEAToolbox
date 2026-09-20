function ok = i_openoutputfile(outfile)
%I_OPENOUTPUTFILE Open a generated output file when possible.
%
%   pkg.i_openoutputfile(outfile)
%   ok = pkg.i_openoutputfile(outfile)
%
%   OK is false when the file is missing or the launch threw. It is not a
%   report that the file is now on screen: on Windows neither RPTVIEW nor
%   WINOPEN raises when the shell declines to open a file, so OK comes
%   back true whether PowerPoint appeared or nothing happened at all.
%
%   So do not rely on this to tell the user their output arrived. Write
%   to a path they chose, and the question stops mattering -- a deck that
%   fails to launch is still a deck they can find. The warning below is
%   not enough on its own either: a caller may have warnings off.

if nargout > 0, ok = false; end

if nargin < 1 || isempty(outfile) || ~isfile(outfile)
    return;
end

[~, ~, ext] = fileparts(outfile);
ext = lower(ext);

try
    switch ext
        case '.docx'
            rptview(outfile, 'docx');
        case '.pptx'
            rptview(outfile);
        otherwise
            if ispc
                winopen(outfile);
            else
                open(outfile);
            end
    end
    if nargout > 0, ok = true; end
catch ME
    warning('Could not open output file "%s": %s', outfile, ME.message);
end

end
