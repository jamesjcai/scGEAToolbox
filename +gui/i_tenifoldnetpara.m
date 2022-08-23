function [nsubsmpl,csubsmpl,savegrn]=i_tenifoldnetpara

%   addOptional(p,'nsubsmpl',10,@(x) fix(x)==x & x>0);
%   addOptional(p,'csubsmpl',500,@(x) fix(x)==x & x>0);
%   addOptional(p,'savegrn',false,@islogical);
nsubsmpl=[];
csubsmpl=[];
savegrn=[];

definput = {'10','500'};
prompt = {'Number of subsamples (nsubsmpl=[10..50]:',...
          'Number of cells per subsample (csubsmpl=[200..5000]):'};
dlgtitle = 'scTenifoldNet Settings';
dims = [1 55];
answer = inputdlg(prompt,dlgtitle,dims,definput);

if isempty(answer), return; end
try
    nsubsmpl=str2double(answer{1});
    csubsmpl=str2double(answer{2});
    assert((nsubsmpl>=10) && (nsubsmpl<=50));
    assert((csubsmpl>=200) && (nsubsmpl<=5000));
catch
    errordlg('Invalid parameter value(s).');
    return;
end

answer=questdlg('Save constructed networks (to current folder)?');
switch answer
    case 'Yes'
        savegrn=true;
    case 'No'
        savegrn=false;
    otherwise
        return;
end

