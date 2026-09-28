function resolution = i_askresolution(parentfig, defaultResolution, prompt)
%I_ASKRESOLUTION Ask for a Louvain resolution, a positive number.
%
%   resolution = gui.i_askresolution(parentfig) asks with Seurat's
%   FindClusters default, 0.8, filled in. DEFAULTRESOLUTION and PROMPT
%   replace the default and the question. RESOLUTION is [] when the user
%   cancels or types something that is not a positive number, after saying
%   which.
%
%   See also GUI.CALLBACK_RECLUSTERCELLS, GUI.CALLBACK_SUBCLUSTERGROUP.

if nargin < 1, parentfig = []; end
if nargin < 2 || isempty(defaultResolution), defaultResolution = 0.8; end
if nargin < 3 || isempty(prompt)
    prompt = 'Louvain resolution (> 0; Seurat default 0.8):';
end

resolution = [];
definput = string(defaultResolution);
if gui.i_isuifig(parentfig)
    answer = gui.i_inputdlg(prompt, definput, parentfig);
else
    answer = inputdlg(prompt, '', [1, 50], {char(definput)});
end
if isempty(answer), return; end
value = str2double(answer{1});
if ~isfinite(value) || value <= 0
    gui.myErrordlg(parentfig, ...
        'The resolution must be a positive number, such as 0.8.', '');
    return;
end
resolution = value;
end
