function [nsubsmpl, csubsmpl, savegrn] = i_tenifoldnetpara(parentfig)

if nargin<1, parentfig = []; end
%   addOptional(p,'nsubsmpl',10,@(x) fix(x)==x & x>0);
%   addOptional(p,'csubsmpl',500,@(x) fix(x)==x & x>0);
%   addOptional(p,'savegrn',false,@islogical);
nsubsmpl = [];
csubsmpl = [];
savegrn = [];

definput = {'10', '500'};
prompt = {'Number of subsamples (nsubsmpl=[10..50]):', ...
    'Number of cells per subsample (csubsmpl=[200..5000]):'};
dlgtitle = 'scTenifoldNet Settings';
dims = [1, 50];
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer), return; end
try
    nsubsmpl = str2double(answer{1});
    csubsmpl = str2double(answer{2});
    assert(isfinite(nsubsmpl) & nsubsmpl==floor(nsubsmpl));
    assert(isfinite(csubsmpl) & csubsmpl==floor(csubsmpl));    
    assert((nsubsmpl >= 10) && (nsubsmpl <= 50));
    assert((csubsmpl >= 200) && (nsubsmpl <= 5000));
catch
    gui.myErrordlg(parentfig, 'Invalid parameter value(s).');
    return;
end

answer = gui.myQuestdlg(parentfig, 'Save constructed networks (to current folder)?');
switch answer
    case 'Yes'
        savegrn = true;
    case 'No'
        savegrn = false;
    otherwise
        return;
end
