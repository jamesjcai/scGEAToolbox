function [cs,ctselected]=callback_CellTypeMarkerScores(src,~,sce)
    if nargin<3
        FigureHandle=src.Parent.Parent;
        sce=guidata(FigureHandle);
    end
cs=[];
ctselected=[];
species=questdlg('Which species?','Select Species','Human','Mouse','Mouse');
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
% idx=find(matches(ctlist,ctselected));
idx=matches(ctlist,ctselected);
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

    if nargout==0
        f0=figure('Visible',false);
        gui.i_stemscatter(sce.s,cs);
        zlabel('Score Value')
        title(ctselected)

tb = uitoolbar(f0);
pkg.i_addbutton2fig(tb,'off',@i_saveCrossTable,"export.gif",'Save cross-table');
pkg.i_addbutton2fig(tb,'off',{@gui.i_savemainfig,3},"powerpoint.gif",'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb,'on',@gui.i_pickcolormap,'plotpicker-compass.gif','Pick new color map...');
pkg.i_addbutton2fig(tb,'on',@gui.i_invertcolor,'plotpicker-comet.gif','Invert colors');
pkg.i_addbutton2fig(tb,'on',@i_geneheatmapx,'plotpicker-cometx.gif','Heatmap');
movegui(f0,'center');
set(f0,'Visible',true);
    end

    function i_saveCrossTable(~,~)
        gui.i_exporttable(cs,false,'CellScore');
    end
    function i_geneheatmapx(~,~)
        gui.i_geneheatmap(sce,sce.c_cell_type_tx,posg);
    end
    
end
