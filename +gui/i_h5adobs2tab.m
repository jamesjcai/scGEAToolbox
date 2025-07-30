function [T] = i_h5adobs2tab(filenm, parentfig)
if nargin < 2, parentfig = []; end
T = [];

hinfo = h5info(filenm);
idx = find(string({hinfo.Groups.Name})=="/obs");
if isempty(idx), return; end

% dnames = {hinfo.Groups(idx).Datasets.Name}
dnames = pkg.i_extractfield(hinfo.Groups(idx).Datasets, 'Name');
if isempty(dnames), return; end

listitems = string(dnames);

if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, listitems, ...
        'Select /obs/Variable:', listitems, true);
else                    
    [indx, tf] = listdlg('PromptString', 'Select /obs/Variable', ...
        'SelectionMode', 'multiple', 'ListString', ...
        listitems, 'ListSize', [260, 300]);
end
if tf ~= 1, return; end

T = table();
for k = 1:length(indx)
    fieldnm = dnames{indx(k)};
    % fieldnm = ['/obs/' dnames{indx(k)}];
    T.(fieldnm) = h5read(filenm, ...
        ['/obs/' fieldnm]);
end
% gui.TableViewerApp(T, parentfig);
