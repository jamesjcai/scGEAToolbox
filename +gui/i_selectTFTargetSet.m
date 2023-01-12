function [glist,setname]=i_selectTFTargetSet(species)
if nargin<1, species='human'; end
glist=[];
setname=[];
pw1=fileparts(mfilename('fullpath'));
switch lower(species)
    case {'human','hs'}        
        f='dorothea_hs.mat';
    case {'mouse','mm'}
        f='dorothea_mm.mat';
end
pth=fullfile(pw1,'..','resources','DoRothEA_TF_Target_DB',f);
load(pth,'T');
[t]=crosstab(T.tf,T.target);   % TF-by-target regulagory relationship matrix
[gid,gnlist]=grp2idx(T.target);
[tid,tflist]=grp2idx(T.tf);


[indx1,tf1] = listdlg('PromptString',...
    {'Select TF gene:'},...
     'SelectionMode','single','ListString',tflist);

if tf1~=1, return; end
    glist=string(gnlist(gid(t(tid(indx1),:)>0)));
    setname=tflist{indx1};
end



