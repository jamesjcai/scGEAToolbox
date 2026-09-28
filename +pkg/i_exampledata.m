function f = i_exampledata(name)
% I_EXAMPLEDATA Local path to an example data file, downloading if needed
%
%   F = PKG.I_EXAMPLEDATA(NAME) returns the full path of NAME in the
%   toolbox's example_data folder. The packaged toolbox ships only
%   new_example_sce.mat, so any other file is downloaded on first use from
%   the public scGEAToolbox repository. The copy is saved in example_data
%   when that folder is writable, otherwise under tempdir, and reused on
%   later calls.
%
%   Example:
%     [X, g] = sc_readfile(pkg.i_exampledata("GSM3204304_P_P_Expr.csv"));

arguments
    name (1,1) string
end

pw1 = fileparts(mfilename('fullpath'));
f = fullfile(pw1, '..', 'example_data', name);
if isfile(f), return; end
fcache = fullfile(tempdir, 'scGEAToolbox', 'example_data', name);
if isfile(fcache), f = fcache; return; end

url = "https://github.com/jamesjcai/scGEAToolbox/raw/refs/heads/main/example_data/" + name;
options = weboptions('Timeout', 60);
fprintf('Downloading example data %s...\n', name);
try
    websave(f, url, options);
catch
    % example_data is read-only in some installs (MATLAB Online, shared
    % toolbox folders); fall back to a per-user cache.
    if isfile(f), delete(f); end
    if ~isfolder(fileparts(fcache)), mkdir(fileparts(fcache)); end
    f = fcache;
    websave(f, url, options);
end
end
