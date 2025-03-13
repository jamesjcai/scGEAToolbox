function [glist, setname] = ...
    i_selectTFTargetSet(species, parentfig)
if nargin < 2, parentfig=[]; end
if nargin < 1, species = 'human'; end
glist = [];
setname = [];
pw1 = fileparts(mfilename('fullpath'));
switch lower(species)
    case {'human', 'hs'}
        f = 'dorothea_hs.mat';
    case {'mouse', 'mm'}
        f = 'dorothea_mm.mat';
end
pth = fullfile(pw1, '..', 'resources', 'DoRothEA_TF_Target_DB', f);
fprintf('\nReading ... %s.\n', pth);
load(pth, 'T');
T = T(T.mor > 0, :); % only consider positive regulation
%size(T)
[t] = crosstab(T.tf, T.target); % TF-by-target regulagory relationship matrix
%size(t)
%assignin('base','t1',t);
[~, gnlist] = grp2idx(T.target);
[~, tflist] = grp2idx(T.tf);

if gui.i_isuifig(parentfig)
    [indx1, tf1] = gui.myListdlg(parentfig, tflist, 'Select TF gene:');
else
    [indx1, tf1] = listdlg('PromptString', ...
        {'Select TF gene:'}, ...
        'SelectionMode', 'single', 'ListString', tflist, 'ListSize', [220, 300]);
end
if tf1 ~= 1, return; end
%indx1
sum((t(indx1, :) ~= 0))
glist = string(gnlist(t(indx1, :) ~= 0));
% glist=string(gnlist(gid(t(tid(indx1),:)>0)));
setname = tflist{indx1};
end
