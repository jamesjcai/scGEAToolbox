function callback_CellTypeMarkerScores(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);


    species=questdlg('Which species?','Select Species','Mouse','Human','Mouse');
switch lower(species)
    case 'human'
        stag='hs';
    case 'mouse'
        stag='mm';
    otherwise
        return;
end

oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'..','+run','thirdparty','alona_panglaodb');
cd(pth);


markerfile=sprintf('marker_%s.mat',stag);
if exist(markerfile,'file')
    load(markerfile,'Tw','Tm');
else
    Tw=readtable(sprintf('markerweight_%s.txt',stag));
    Tm=readtable(sprintf('markerlist_%s.txt',stag),...
        'ReadVariableNames',false,'Delimiter','\t');
    % save(markerfile,'Tw','Tm');
end
cd(oldpth);

ctlist=string(Tm.Var1);
listitems=sort(ctlist);
[indx,tf] = listdlg('PromptString',...
    {'Select Class','',''},...
     'SelectionMode','single','ListString',listitems,'ListSize',[220,300]);
if ~tf==1, return; end
ctselected=listitems(indx);
idx=find(matches(ctlist,ctselected));
ctmarkers=Tm.Var2{idx};

posg=string(strsplit(ctmarkers,','));
posg(strlength(posg)==0)=[];
cs=sc_cellscore(sce.X,sce.g,posg);

posg=sort(posg);
fprintf('\n=============\n%s\n-------------\n',ctselected);
for k=1:length(posg)
    fprintf('%s\n',posg(k));
end
fprintf('=============\n');
figure;
gui.i_stemscatter(sce.s,cs);
zlabel('Score Value')
title(ctselected)
end
