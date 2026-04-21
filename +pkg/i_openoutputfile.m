function i_openoutputfile(outfile)
%I_OPENOUTPUTFILE Open a generated output file when possible.

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
catch ME
    warning('Could not open output file "%s": %s', outfile, ME.message);
end

end
