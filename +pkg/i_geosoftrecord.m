function [lines, url] = i_geosoftrecord(acc)
% I_GEOSOFTRECORD  Download a GEO record as SOFT-format text lines.
%
%   [lines, url] = pkg.i_geosoftrecord(acc)
%
%   Inputs:
%     acc   - GEO accession ('GSM...', 'GSE...', 'GPL...'), or a full
%             acc.cgi URL carrying acc=<accession>
%
%   Outputs:
%     lines - string array holding one SOFT line per element: a leading
%             '^ENTITY = ACC' line followed by '!Entity_attribute = value'
%     url   - the URL the record was read from
%
%   NCBI serves the rendered acc.cgi HTML page behind a CAPTCHA for
%   automated clients, so scraping that page yields an empty record with no
%   HTTP error to give the failure away. The SOFT text view is not gated,
%   and it is already machine-readable.
%
%   See also pkg.i_geosoftfield, sc_readgeoaccess

    acc = string(acc);
    if startsWith(acc, "http", 'IgnoreCase', true)
        tok = regexp(acc, 'acc=([^&]+)', 'tokens', 'once');
        if isempty(tok)
            error('pkg:i_geosoftrecord:NoAccession', ...
                'No GEO accession found in %s.', acc);
        end
        acc = string(tok(1));
    end
    acc = upper(strtrim(acc));

    url = sprintf(['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi' ...
        '?acc=%s&targ=self&form=text&view=brief'], acc);

    % NCBI answers with a 502 or simply stalls often enough that webread's
    % 5-second default turns routine flakiness into a hard failure. Give it
    % room, and take one more run at it before giving up.
    opts = weboptions('Timeout', 30);
    try
        raw = webread(url, opts);
    catch
        pause(2);
        raw = webread(url, opts);
    end
    lines = string(splitlines(string(raw)));

    % A SOFT record always opens with ^SAMPLE / ^SERIES / ^PLATFORM. Anything
    % else means GEO answered with the HTML page instead of the record.
    if ~any(startsWith(lines, "^"))
        error('pkg:i_geosoftrecord:NoRecord', ...
            ['GEO did not return a record for %s. Check the accession, ' ...
            'or try again later if NCBI is throttling requests.'], acc);
    end
end
