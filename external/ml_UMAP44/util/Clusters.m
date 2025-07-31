classdef Clusters < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(SetAccess=private)
        ids;
        clues;
        sizes;
        orderIdxs;
        fcsIdxs;
        N;
    end
    properties
        colors;
        selColors;
    end
    
    methods
        function this=Clusters(clues, cluColorFactory)
            ids=unique(clues);
            zeroIdx=find(ids==0,1);
            if ~isempty(zeroIdx)
                ids(zeroIdx)=[];
            end
            sizes=MatBasics.HistCounts(clues, ids);
            orderIdxs=MatBasics.GetSizeRankings(sizes);
            N=length(ids);
            if nargin>1
                if isempty(zeroIdx)
                    colors=zeros(N, 3);
                    selColors=zeros(N, 3);
                    ii=0;
                else
                    colors=zeros(N+1, 3);
                    selColors=zeros(N+1, 3);
                    colors(1,:)=[0 0 0];
                    selColors(1,:)=[0 0 0];
                    ii=1;
                end
                for i=1:N
                    colors(i+ii, :)=cluColorFactory.get(orderIdxs(i), false);
                    selColors(i+ii, :)=cluColorFactory.get(orderIdxs(i), true);
                end
                this.colors=colors;
                this.selColors=selColors;
            end
            fcsIdxs=cell(1, N);
            for i=1:N
                id=ids(i);
                fcsIdxs{i}=find(clues==id);
            end
            this.fcsIdxs=fcsIdxs;
            this.ids=ids;
            this.clues=clues;
            this.sizes=sizes;
            this.orderIdxs=orderIdxs;
            this.N=N;
        end
        
        function f=fMeasure(this, thisI, that, thatI)
            same=intersect(this.fcsIdxs{thisI}, that.fcsIdxs{thatI});
            f=Clusters.F_measure(length(same), this.sizes(thisI), ...
                that.sizes(thatI));
        end
    end
    
    methods(Static)
        function [f1Score, precision, recall]=F_measure(...
                sizeTruePos, sizeTrainingSet, sizeTestSet)
            precision=sizeTruePos/sizeTestSet;
            recall=sizeTruePos/sizeTrainingSet;
            f1Score=(2*precision*recall)/(precision+recall);
            if isnan(f1Score)
                f1Score=0;
            end
        end
       
         function f=F_measures(nsTruePos, sizeTrainingSet, sizesTestSet)
            precision=nsTruePos./sizesTestSet;
            recall=nsTruePos./sizeTrainingSet;
            f=(2*precision.*recall)./(precision+recall);
            if any(isnan(f))
                f(isnan(f))=0;
            end
         end
         
         function p=Precision(nTruePos, sizesTestSet)
            p=nTruePos./sizesTestSet;
            if any(isnan(p))
                p(isnan(p))=0;
            end
         end
         
         function r=Recall(nTruePos, sizeTrainingSet)
            r=nTruePos./sizeTrainingSet;
            if any(isnan(r))
                r(isnan(r))=0;
            end
         end

         function c=Concordance(nTruePos, sizeTrainingSet, sizeTestSet)
            if sizeTrainingSet == 0 && sizeTestSet == 0
                c = NaN;
            else
                c = nTruePos/(sizeTrainingSet + sizeTestSet - nTruePos);
            end
         end

         function c=Concordances(nsTruePos, sizeTrainingSet, sizesTestSet)
            c=nsTruePos./(sizeTrainingSet + sizesTestSet - nsTruePos);
            c(isnan(c))=0;
         end
         
         function [clrs, selClrs]=GrayColors(clustSizes)
             numClusts=length(clustSizes);
             clrs=zeros(numClusts+1,3);
             selClrs=zeros(numClusts+2,3);
             [~,I]=sort(clustSizes, 'descend');
             grayGap=.8/numClusts;
             for idx_=2:numClusts
                 cluIdx=I(idx_);
                 %cluIdx=idx_;
                 clr=.04+(idx_*grayGap);
                 rank=mod(idx_, 3);
                 if rank==0
                     greenFactor=1;
                     blueFactor=1.1;
                 elseif rank==1
                     greenFactor=1.15;
                     blueFactor=1;
                 else
                     blueFactor=.98;
                     greenFactor=1.03;
                 end
                 clrs(cluIdx,:)=[clr greenFactor*clr clr ];
                 selClrs(cluIdx,:)=.9*clrs(cluIdx,:);
                 %clrs(cluIdx,:)=[clr clr*greenFactor blueFactor*clr ];
             end
         end
         
         function F=MeasureF(this, that)
            nThis=length(this.clues);
            F=0;
            for i=1:this.N
                a=zeros(1, that.N);
                for j=1:that.N
                    a(j)=this.fMeasure(i, that, j);
                end
                F=F+(this.sizes(i)/nThis*max(a));
            end
         end

         function [f1s, sizeIntersections, sizeTrainingSets, ...
                sizeTestSets]=F1_scores(data, ...
                trLabels, teLabels, trainingLabels)
            if isempty(data)
                return;
            end
            if nargin<4
                trainingLabels=unique(trLabels);
            end
            if isnumeric(trainingLabels)
                trainingLabels=trainingLabels(trainingLabels~=0);
            end
            R=size(data,1);
            if R ~= length(trLabels) ...
                    || R ~= length(teLabels)
                f1s=[];
                sizeIntersections=[]; 
                sizeTrainingSets=[]; ...
                sizeTestSets=[];
                warning(['# of classification IDs must equal # ' ...
                    'of data points']);
                return;
            end
            try
                N=length(trainingLabels);
                f1s=zeros(N, 1);
                sizeTrainingSets=zeros(N, 1);
                sizeIntersections=zeros(N, 1);
                sizeTestSets=zeros(N, 1);
                for i=1:N
                    if iscell(trLabels)
                        label=trainingLabels{i};
                        trChoices=strcmp(trLabels, label);
                        teChoices=strcmp(teLabels, label);
                    else
                        label=trainingLabels(i);
                        trChoices=trLabels==label;
                        teChoices=teLabels==label;
                    end
                    sizeIntersections(i)=sum(trChoices&teChoices);
                    sizeTrainingSets(i)=sum(trChoices);
                    sizeTestSets(i)=sum(teChoices);
                    if sizeTrainingSets(i) == 0 || sizeTestSets(i) == 0
                        f1s(i) = 0;
                    else
                        precision=sizeIntersections(i)/sizeTestSets(i);
                        recall=sizeIntersections(i)/sizeTrainingSets(i);
                        f1s(i)=(2*precision.*recall)/(precision+recall);
                        if isnan(f1s(i))
                            f1s(i)=0;
                        end
                    end
                end
            catch ex
                ex.getReport
            end
         end

         function similarities=Similarities(trData, trLabels, ...
                 teData, teLabels, contextName, trainingLabels, ...
                 compareTo, debug)
             if nargin<8
                 debug=true;
                 if nargin<7
                     compareTo=[];
                     if nargin<6
                         trainingLabels=unique(trainingLabels);
                         if nargin<5
                             contextName='';
                         end
                     end
                 end
             end
             if ~check('training set', trData, trLabels) ||...
                 ~check('test set', teData, teLabels) 
                 similarities=[];
                 return;
             end
             if debug
                 fprintf(['Creating probability bins ' ...
                     'for %dx%d and %dx%d\n'], size(trData, 1),...
                     size(trData, 2), size(teData, 1),...
                     size(teData, 2));
             end
             ab=SuhProbabilityBins.Bins(trData, teData);
             isCharLabels=iscell(trainingLabels);
             if ~isCharLabels
                 %eliminate background 0
                 trainingLabels=trainingLabels(trainingLabels~=0);
             end
             N=length(trainingLabels);
             similarities=zeros(N,1);
             for i=1:N
                 if isCharLabels
                     label=trainingLabels{i};
                 else
                     label=trainingLabels(i);
                 end
                 trRows=strcmp(trLabels, label);
                 teRows=strcmp(teLabels, label);
                 %assign to scalar 1st for ease of debug watching
                 similarity=1-ab.distance(trRows, teRows, false);
                 if similarity<0
                     similarity=0; %-99 is used for legacy reasons
                 end
                 if debug
                     if ~isempty(compareTo)
                         dif=compareTo(i)-similarity;
                         if abs(dif)>.01
                             fprintf('***');
                         end
                     end
                     if isCharLabels
                         fprintf(['%s %s similarity is %s between ' ...
                             '%d X %d points\n'],contextName, label, ...
                             String.encodePercent(similarity), sum(trRows), ...
                             sum(teRows));
                     else
                         fprintf(['%s %d similarity is %s between ' ...
                             '%d X %d points\n'],contextName, label, ...
                             String.encodePercent(similarity), sum(trRows), ...
                             sum(teRows));

                     end
                 end
                 similarities(i)=similarity;
             end

             function ok=check(name, data, labels)
                 nD=size(data, 1);
                 nL=length(labels);
                 if nD ~= nL
                     warning(['NOTHING to do %s data rows=%d ' ...
                         'and # of labels =%d:'], name, nD, nL);
                     ok=false;
                 else
                     ok=true;
                 end
             end
         end

         function [iji, ijiSizeIntersection, ijiSizeTrainingSets, ...
                ijiSizeTestSets, ji, jiSizeIntersection, jiSizeTrainingSets, ...
                jiSizeTestSets]=InlierJaccardIndex(data, ...
                trLabels, teLabels, trainingLabels, coreCutoff)
            if isempty(data)
                return;
            end
            if nargin<5
                coreCutoff=.8;
                if nargin<4
                    trainingLabels=unique(trLabels);
                end
            end
            if isnumeric(trainingLabels)
                trainingLabels=trainingLabels(trainingLabels~=0);
            end
            if coreCutoff>=1
                coreCutoff = 0.8;
            end
            R=size(data,1);
            if R ~= length(trLabels) ...
                    || R ~= length(teLabels)
                iji=[];
                ji=[];
                ijiSizeIntersection=[]; 
                ijiSizeTrainingSets=[]; ...
                ijiSizeTestSets=[];
                warning(['# of classification IDs must equal # ' ...
                    'of data points']);
                return;
            end
            try
                N=length(trainingLabels);
                iji=zeros(N, 1);
                ji=zeros(N, 1);
                ijiSizeTrainingSets=zeros(N, 1);
                ijiSizeIntersection=zeros(N, 1);
                ijiSizeTestSets=zeros(N, 1);
                jiSizeTrainingSets=zeros(N, 1);
                jiSizeIntersection=zeros(N, 1);
                jiSizeTestSets=zeros(N, 1);
                for i=1:N
                    %training set choices of inliers
                    if iscell(trLabels)
                        label=trainingLabels{i};
                        trChoices=strcmp(trLabels, label);
                        teChoices=strcmp(teLabels, label);
                    else
                        label=trainingLabels(i);
                        trChoices=trLabels==label;
                        teChoices=teLabels==label;
                    end
                    ji(i)=MatBasics.Concordance(trChoices, teChoices);
                    coreUnionChoices=trChoices|teChoices;
                    jiSizeIntersection(i)=sum(trChoices&teChoices);
                    jiSizeTrainingSets(i)=sum(trChoices);
                    jiSizeTestSets(i)=sum(teChoices);
                    [core,ex]=MatBasics.MahalanobisCore( ...
                        data(coreUnionChoices,:),coreCutoff);
                    if ~isempty(ex)
                        warning('"%s" training set %d events had exception\n%s', ...
                            label, sum(trChoices), ex.message)
                    end
                    coreUnionChoices(coreUnionChoices) = core;
                    sizeRecalled=sum(trChoices&teChoices&coreUnionChoices);
                    sizeTrainingSet=sum(trChoices&coreUnionChoices);
                    sizeTestSet=sum(teChoices&coreUnionChoices);
                    ijiSizeIntersection(i)=sizeRecalled;
                    ijiSizeTrainingSets(i)=sum(trChoices & coreUnionChoices);
                    ijiSizeTestSets(i)=sum(teChoices & coreUnionChoices);
                    if ijiSizeTrainingSets(i) == 0 || ijiSizeTestSets(i) == 0
                        iji(i) = 0;
                    else
                        iji(i)=Clusters.Concordance(sizeRecalled, sizeTrainingSet, sizeTestSet);
                    end
                end
            catch ex
                ex.getReport
            end
         end

         %cube is test cases X cell types X measurements per cell type
         function similarities=SimilarityCube(cube, fcn)
             if nargin<2
                 fcn=@(x)median(x);
             end
             similarities=feval(fcn, (cube(:,:,1)))';
         end

         function tbl=SimilarityTable(cubes, cubeNames, header1, column1, fcn)
             if nargin<5
                 fcn=@(x)median(x);
             end
             N=length(cubes);
             varArg=cell(1, 1+N);
             varArg{1}=column1;
             varNames=cell(1, 1+N);
             varNames{1}=[header1 ' ' char(fcn)];
             for i=1:N
                 cubeName=cubeNames{i};
                 idx=2+(i-1);
                 varArg{idx}=Clusters.SimilarityCube(cubes{i}, fcn);
                 varNames{idx}=[cubeName ' Similarity'];
             end
             varArg{end+1}='VariableNames';
             varArg{end+1}=varNames;
             tbl=table(varArg{:});
         end

         %cube is test cases X cell types X measurements per cell type
         function [ics, intersections, unions]=IcCube(cube, fcn)
             if nargin<2
                 fcn=@(x)median(x);
             end
             ics=feval(fcn, (cube(:,:,1)))';
             intersections=round(feval(fcn, (cube(:,:,2))))';
             unions=round(feval(fcn, (cube(:,:,3)))'...
                 +feval(fcn, (cube(:,:,4)))'-intersections);
         end

         function tbl=IcTable(cubes, cubeNames, header1, column1, fcn)
             if nargin<5
                 fcn=@(x)median(x);
             end
             N=length(cubes);
             varArg=cell(1, 1+(N*3));
             varArg{1}=column1;
             varNames=cell(1, 1+(N*3));
             varNames{1}=[header1 ' ' char(fcn)];
             for i=1:N
                 cubeName=cubeNames{i};
                 idx=2+(i-1);
                 [varArg{idx}, varArg{idx+N}, varArg{idx+(N*2)}]...
                     =Clusters.IcCube(cubes{i}, fcn);
                 varNames{idx}=[cubeName ' IJI'];
                 varNames{idx+N}=[cubeName ' intersection'];
                 varNames{idx+(N*2)}=[cubeName ' union'];
             end
             varArg{end+1}='VariableNames';
             varArg{end+1}=varNames;
             tbl=table(varArg{:});
         end

         function tbl=JiTable(cubes, cubeNames, header1, column1, fcn)
             if nargin<5
                 fcn=@(x)median(x);
             end
             N=length(cubes);
             varArg=cell(1, 1+(N*3));
             varArg{1}=column1;
             varNames=cell(1, 1+(N*3));
             varNames{1}=[header1 ' ' char(fcn)];
             for i=1:N
                 cubeName=cubeNames{i};
                 idx=2+(i-1);
                 [varArg{idx}, varArg{idx+N}, varArg{idx+(N*2)}]...
                     =Clusters.IcCube(cubes{i}, fcn);
                 varNames{idx}=[cubeName ' JI'];
                 varNames{idx+N}=[cubeName ' intersection'];
                 varNames{idx+(N*2)}=[cubeName ' union'];
             end
             varArg{end+1}='VariableNames';
             varArg{end+1}=varNames;
             tbl=table(varArg{:});
         end
         %cube is test cases X cell types X measurements per cell type
         function [f1s, intersections, trSizes, teSizes]=F1Cube(cube, fcn)
             if nargin<2
                 fcn=@(x)median(x);
             end
             f1s=feval(fcn, (cube(:,:,1)))';
             intersections=round(feval(fcn, (cube(:,:,2))))';
             trSizes=round(feval(fcn, (cube(:,:,3))))';
             teSizes=round(feval(fcn, (cube(:,:,4))))';
         end

         function tbl=F1Table(cubes, cubeNames, header1, column1, fcn)
             if nargin<5
                 fcn=@(x)median(x);
             end
             N=length(cubes);
             varArg=cell(1, 1+(N*4));
             varArg{1}=column1;
             varNames=cell(1, 1+(N*4));
             varNames{1}=[header1 ' ' char(fcn)];
             for i=1:N
                 cubeName=cubeNames{i};
                 idx=2+(i-1);
                 [varArg{idx}, varArg{idx+N}, varArg{idx+(N*2)}, ...
                     varArg{idx+(N*3)}]=Clusters.F1Cube(cubes{i}, fcn);
                 varNames{idx}=[cubeName ' F1'];
                 varNames{idx+N}=[cubeName ' intersection'];
                 varNames{idx+(N*2)}=[cubeName ' # trainers'];
                 varNames{idx+(N*3)}=[cubeName ' # tests'];
             end
             varArg{end+1}='VariableNames';
             varArg{end+1}=varNames;
             tbl=table(varArg{:});
         end
    end
end