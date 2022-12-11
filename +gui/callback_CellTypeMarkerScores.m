function [cs,ctselected]=callback_CellTypeMarkerScores(src,~,sce)
    if nargin<3
        FigureHandle=src.Parent.Parent;
        sce=guidata(FigureHandle);
    end
cs=[];
ctselected=[];
species=questdlg('Which species?','Select Species','Human','Mouse','Human');
switch lower(species)
    case 'human'
        stag='hs';
    case 'mouse'
        stag='mm';
    otherwise
        return;
end


actiontype=questdlg('Cell score for PanglaoDB marker genes or MSigDB genes?',...
    '','PanglaoDB Markder Genes',...
    'MSigDB Genes','PanglaoDB Markder Genes');

switch actiontype
    case 'PanglaoDB Markder Genes'
    
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
        {'Select Class'},...
         'SelectionMode','single','ListString',listitems,'ListSize',[220,300]);
        if ~tf==1, return; end
        ctselected=listitems(indx);
        % idx=find(matches(ctlist,ctselected));
        idx=matches(ctlist,ctselected);
        ctmarkers=Tm.Var2{idx};
        posg=string(strsplit(ctmarkers,','));
        posg(strlength(posg)==0)=[];

    case 'MSigDB Genes'
        [posg,ctselected]=gui.i_selectMSigDBGeneSet(stag);
        if isempty(posg) || isempty(ctselected)
            return;
        end
end

    
    % matches(sce.g, posg,'IgnoreCase',true);

    [cs]=sc_cellscore(sce.X,sce.g,posg);
    
    posg=sort(posg);
    fprintf('\n=============\n%s\n-------------\n',ctselected);
    for k=1:length(posg)
        fprintf('%s\n',posg(k));
    end
    fprintf('=============\n');
    if nargout==0
        gui.i_stemscatterfig(sce,cs,posg,ctselected);
        % a=inputdlg('Gene Set Info:','Gene Viewer',[10 50],{char(sce.metadata)});
    end



%     function i_saveCrossTable(~,~)
%         gui.i_exporttable(cs,false,'CellScore');
%     end
%     function i_geneheatmapx(~,~)
%         gui.i_geneheatmap(sce,sce.c_cell_type_tx,posg);
%     end
    
end
