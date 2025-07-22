function [qf, qfTree, fig]=run_QfTree(trainingSet, trainingIds, ttl, varargin)
%   Class for Hi-D matching with merging of data subsets using QF match or
%   F-measure or both.

%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause
%
%   QF Tree (phenogram) is described in 
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/

ml=BasicMap.Global.getNumeric(QfHiDM.PROP_MERGE_LIMIT, 1);
p=parseArguments();
parse(p,varargin{:});
[args, hadPu]=checkArgs(p.Results);
if args.askTitle
    nLines=length(ttl)+1;
    labels=cell(1, nLines);
    ttl{end+1}='';
    labels{1}='N subsets in...';
    for ii=2:nLines
        labels{ii}=['Line #' num2str(ii)];
    end
    if ~isempty(args.locate_fig) && iscell(args.locate_fig)
        jw=Gui.JWindow(args.locate_fig{1});
    else
        jw=[];
    end
    inputs=inputsDlg('Enter title for phenogram', 'Confirm...', ...
        labels, ttl, 'center', false, 20, 2, 1,  [], [], [], jw, ...
        false, [], [], [], false);
    if isempty(inputs)
        qf=[]; qfTree=[]; fig=[];
        if ~hadPu
            args.pu.close;
        end
        return;
    end
    ttl={};
    for ii=1:nLines
        if ~isempty(inputs{ii})
            ttl{end+1}=inputs{ii};
        end
    end
end
qf=QfHiDM(trainingSet, args.trainingSetComp, trainingIds, ...
    trainingSet, args.trainingSetComp, trainingIds, args.bins, args.binStrategy);
qf.tNames=args.trainingNames;
qf.sNames=args.trainingNames;
qf.matchStrategy=args.matchStrategy;
qf.mergeStrategy=args.mergeStrategy;
qf.maxDeviantParameters=0;
qf.computeTree(args.pu, qf.tNames);
qfTree=QfTree(qf, ttl, BasicMap.Global, '', args.colors, ...
    args.edgeColors, args.lineWidths, qf.tNames);
fig=qfTree.fig;
if ~isempty(args.locate_fig)
    if islogical(args.locate_fig) && args.locate_fig
        set(fig, 'visible', 'on');
    elseif iscell(args.locate_fig)
        SuhWindow.Follow(fig, args.locate_fig);
        SuhWindow.SetFigVisible(fig);
    end
end
if ~hadPu
    args.pu.close;
end
BasicMap.Global.setNumeric(QfHiDM.PROP_MERGE_LIMIT, num2str(ml));
Gui.AddSvgToToolBar(fig);

    function [args, hadPu]=checkArgs(args)
        hadPu=~isempty(args.pu);
        [trainers, trainingIds]=checkSet(trainingSet, trainingIds, args.trainingNames);
        if args.log10
            trainingSet=QfHiDM.Log10(trainingSet);
        end
        if isempty(args.trainingSetComp)
            args.trainingSetComp=trainingSet;
        end
        if isempty(args.pu)
            args.pu=PopUp(['Building tree for ' num2str(length(trainers)) ...
                ' subsets '], 'north', 'Running QF Tree', false, true);
            args.pu.setTimeSpentTic(tic);
        elseif isequal('none', args.pu)
            args.pu=[];
        end
        BasicMap.Global.setNumeric(QfHiDM.PROP_MERGE_LIMIT, ...
            num2str(args.mergeLimit));
        nNames=length(args.trainingNames);
        if isempty(args.colors)
            args.colors=zeros(nNames, 3);
            for i=1:nNames
                args.colors(i,:)=Gui.HslColor(i, nNames);
            end
        end
        if isempty(args.edgeColors)
            args.edgeColors=args.colors;
        end
        if isempty(args.lineWidths)
            args.lineWidths=zeros(1, nNames)+1;
        end
    end

    function [u, ids]=checkSet(dataSet, ids, names)
        R1=size(dataSet,1);
        [R2, C2]=size(ids);
        assert((R1==R2 && C2==1) || (R1==C2 && R2==1), 'Need same # of data rows and id rows');
        if R1==C2 && R2==1
            ids=ids';
        end
        u=unique(ids(ids ~= 0));
        if ~isempty(names)
            R3=length(names);
            assert(length(u)==R3, 'Need same # of names as non-zero IDs');
        end
    end

    function p=parseArguments(varargin)
        p = inputParser;
        addParameter(p,'bins', 0, @isnumeric);
        addParameter(p,'binStrategy', 0, @isnumeric);
        addParameter(p,'log10', true, @islogical);
        addParameter(p,'askTitle', true, @islogical);
        addParameter(p,'trainingSetComp', [], @isnumeric);
        addParameter(p,'trainingNames', {}, @(x)StringArray.IsOk(x));
        addParameter(p,'pu',[], @(x)isempty(x) || isequal('none', x) || isa(x,'PopUp'));
        addParameter(p,'mergeLimit',3, @(x)validMergeLimit(x));
        addParameter(p,'matchStrategy',1, @(x)validMatchStrategy(x));
        addParameter(p,'mergeStrategy', 1, @(x)validMergeStrategy(x));
        addParameter(p,'colors', [], @(x)isColors(x));
        addParameter(p,'edgeColors', [], @isColors);
        addParameter(p,'lineWidths', [], @isnumeric);
        addParameter(p, 'locate_fig', {}, @(x) islogical(x) ...
            || Args.IsLocateFig(x, 'locate_fig' ));
            
        function ok=validMergeLimit(x)
            ok=x>=1 && x <=12;
        end

        function ok=isColors(x)
            ok=false;
            if isnumeric(x)
                C=size(x, 2);
                if C==3
                    ok=true;
                end
            end
        end
        function ok=validMatchStrategy(x)
            ok=x>=1 && x <=3;
        end
        function ok=validMergeStrategy(x)
            ok=x>=1 && x <=8;
        end

    end

end
