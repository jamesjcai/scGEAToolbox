classdef AdaptiveBins<handle
%   Class for probability binning described at
%   https://www.ncbi.nlm.nih.gov/pubmed/11598946

%   Bioinformatics inventors
%   Original:  Roederer M1, Moore W, Treister A, Hardy RR, Herzenberg LA.
%   Incorporation into QF matching:  
%              Darya Orlova <dyorlova@gmail.com>
%              Guenther Walther <gwalther@stanford.edu>
%              Connor Meehan <connor.gw.meehan@gmail.com>
%
%   Software Developers: Connor Meehan <connor.gw.meehan@gmail.com>
%                        Stephen Meehan <swmeehan@stanford.edu> 
%
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(SetAccess=private)
        means;
        dists;
        teachPtrs;
        studPtrs;
        teachSize;
        studSize;
        maxDist=0;
        binSize;
    end 
    
    properties(Constant)
        MAX_SIZE=25000;
        MIN_SIZE=30;
    end
    
    methods
        function this=AdaptiveBins(teachData, studData, minBinSize, ...
                sizeIsFixedSplit, teachStudCacheFile, studTeachCacheFile)
            if nargin<6
                studTeachCacheFile=[];
                if nargin<5
                    teachStudCacheFile=[];
                    if nargin<4
                        sizeIsFixedSplit=false;
                        if nargin<3
                            minBinSize=[];
                        end
                    end
                end
            end
            if ~isempty(teachStudCacheFile)
                teachStudCacheFile=[teachStudCacheFile 'pb_v3_' num2str(minBinSize) ...
                    '_' num2str(sizeIsFixedSplit) '.mat'];
            end
            if ~isempty(studTeachCacheFile)
                studTeachCacheFile=[studTeachCacheFile 'pb_v3_' num2str(minBinSize) ...
                        '_' num2str(sizeIsFixedSplit) '.mat'];
            end
            if exist(teachStudCacheFile, 'file')
                load(teachStudCacheFile, 'means', 'teachPtrs', 'studPtrs', 'binSize'); 
                this.means=means;
                this.teachPtrs=teachPtrs;
                this.studPtrs=studPtrs;
                this.binSize=binSize;
            elseif exist(studTeachCacheFile, 'file')
                load(studTeachCacheFile, 'means', 'teachPtrs', 'studPtrs', 'binSize'); 
                this.means=means;
                this.teachPtrs=studPtrs;
                this.studPtrs=teachPtrs;
                this.binSize=binSize;
            else
                [this.means, this.teachPtrs, this.studPtrs, testW1, ...
                    testW2, this.binSize]=AdaptiveBins.Create(teachData,...
                    studData, minBinSize, sizeIsFixedSplit, ...
                    teachStudCacheFile);
                assert(round(sum(testW1),11)==1);
                assert(round(sum(testW2),11)==1);
            end
            this.teachSize=size(teachData, 1);
            this.studSize=size(studData, 1);
            if size(this.means, 1)<AdaptiveBins.MAX_SIZE
                this.dists=MatBasics.PDist2Self(this.means);
                this.maxDist=max(max(this.dists));
            end
        end
        
    end
    properties
        uncompressedStudPtrs;
        uncompressedTeachPtrs;
    end
    
    methods
        function compressStudPtrs(this, probability_bins)
            this.uncompressedStudPtrs=this.studPtrs;
            this.studPtrs=probability_bins.fit(this.studPtrs)';
        end
        
        function compressTeachPtrs(this, probability_bins)
            this.uncompressedTeachPtrs=this.teachPtrs;
            this.teachPtrs=probability_bins.fit(this.teachPtrs)';
        end
        
        function uncompress(this)
            if ~isempty(this.uncompressedStudPtrs)
                this.studPtrs=this.uncompressedStudPtrs;
                this.uncompressedStudPtrs=[];
            end
            if ~isempty(this.uncompressedTeachPtrs)
                this.teachPtrs=this.uncompressedTeachPtrs;
                this.uncompressedTeachPtrs=[];
            end
        end
        function [f,p,r]=fMeasure(this, teachChoices, studChoices)
            hPtrs=this.teachPtrs(teachChoices);
            fPtrs=this.studPtrs(studChoices);
            u=unique([hPtrs fPtrs]);
            teachCnts=MatBasics.HistCounts(hPtrs, u);
            studCnts=MatBasics.HistCounts(fPtrs, u);
            [f, p, r]=AdaptiveBins.F_measure(teachCnts,studCnts);
        end
        
        function e=emd(this, teachChoices, studChoices)
            hPtrs=this.teachPtrs(teachChoices);
            fPtrs=this.studPtrs(studChoices);
            u=unique([hPtrs fPtrs]);
            tBinCnts=MatBasics.HistCounts(hPtrs, u);
            sBinCnts=MatBasics.HistCounts(fPtrs, u);
            tSize=this.teachSize;%=sum(teachChoices);
            sSize=this.studSize;%=sum(studChoices);
            tW=tBinCnts/tSize;
            sW=sBinCnts/sSize;
            M=this.means(u, :);
            [~,e]=emd_flow(M, M, tW, sW, @gdf);
        end
        
        function [D, d_max, A_IJ]=distance(this, tChoices, sChoices, weighBySampleSize)
            try
                isMeans=isempty(this.dists);
                [meansOrDists, h, f]=this.weigh(...
                    tChoices, sChoices, weighBySampleSize);
                if isempty(meansOrDists)
                    D = 0;
                    d_max = 0;
                    A_IJ = 0;
                    return
                end
                if nargin<4 || isMeans
                    R=size(meansOrDists, 1);
                    if R>AdaptiveBins.MAX_SIZE
                        if nargin<5 || ~ignoreTooBig
                            D=QfHiDM.ComputeFastHiD(h, f, meansOrDists);
                            if isnan(D)
                                D=QfHiDM.MAX_QF_DISTANCE;%close to max QF distance
                            end
                        else
                            D=QfHiDM.MAX_QF_DISTANCE;
                        end
                        return;
                    end
                    originalMeans=meansOrDists;
                    meansOrDists=MatBasics.PDist2Self(meansOrDists);
                end
                d_max=max(max(meansOrDists));
                A_IJ=1-meansOrDists/d_max;
                v=h-f;
                D=sqrt(v * A_IJ * v');
            catch ex
                MatBasics.SourceLocation(ex)
                if nargin<4 || isMeans
                    D=QfHiDM.ComputeFastHiD(h, f, originalMeans);
                else
                    D=nan;
                end
            end
            if isnan(D)
                D=QfHiDM.MAX_QF_DISTANCE;%close to max QF distance
            end
        end

        function [meansOrDistances, teachWeights, studWeights]=...
                weigh(this, teachChoices, studChoices, weighBySampleSize)
            if weighBySampleSize
                tSize=this.teachSize;%=sum(teachChoices);
                sSize=this.studSize;%=sum(studChoices);
            else
                tSize=sum(teachChoices);
                sSize=sum(studChoices);
            end
            hPtrs=this.teachPtrs(teachChoices);
            fPtrs=this.studPtrs(studChoices);
            u=unique([hPtrs fPtrs]);
            tBinCnts=MatBasics.HistCounts(hPtrs, u);
            sBinCnts=MatBasics.HistCounts(fPtrs, u);
            if AdaptiveBins.DEBUG_LEVEL>0
                f=AdaptiveBins.F_measure(tBinCnts, sBinCnts); %#ok<NASGU> 
            end
            teachWeights=tBinCnts/tSize;
            studWeights=sBinCnts/sSize;
            if QfHiDM.DEBUG_LEVEL>0
                if ~weighBySampleSize
                    assert(round(sum(teachWeights),11)==1);
                    assert(round(sum(studWeights),11)==1);
                end
            end
            if isempty(this.dists)
                meansOrDistances=this.means(u, :);
            else
                meansOrDistances=this.dists(u, u);
                if QfHiDM.DEBUG_LEVEL>0
                    test=MatBasics.PDist2Self(this.means(u, :));
                    if ~isequal(round(test,8), round(meansOrDistances,8))
                        disp('floating point insensitivity at 8th digit?...');
                    end
                end
            end
        end
        
        
    end
    
    methods(Static)
        function [f, precision, recall]=F_measure(teachBinCnts, studBinCnts)
            tSize=sum(teachBinCnts);
            sSize=sum(studBinCnts);
            teachBins=studBinCnts>teachBinCnts&teachBinCnts>0;
            studBins=studBinCnts<=teachBinCnts&studBinCnts>0;
            truePositive=sum(teachBinCnts(teachBins))+sum(studBinCnts(studBins));
            
            % this statement on 5 lines below does nothing much to help
            %truePositive = floor( truePositive + ...
            %    sum((studBinCnts(teachBins)-teachBinCnts(teachBins))...
            %    .* (studBinCnts(teachBins)/sSize)) +...
            %  sum((teachBinCnts(studBins)-studBinCnts(studBins))...
            %    .*  (teachBinCnts(studBins)/tSize)) );
            precision=truePositive/sSize;
            recall=truePositive/tSize;
            f=(2*precision*recall)/(precision+recall);
            if isnan(f)
                f=0;
            end
        end
        
        function [f, p, r]=F_measureOld(teachCnts, studCnts)
            tSize=sum(teachCnts);
            sSize=sum(studCnts);
            cnt=tSize+sSize;
            tRatio=tSize/cnt;
            sRatio=sSize/cnt;
            isOverlap=teachCnts>0 & studCnts>0;
            truePositive=sum(studCnts(isOverlap));
            recalled=sum(teachCnts(isOverlap));
            tRecallRatio=recalled/tSize;
            sRecallRatio=truePositive/sSize;
            ratio=(tRecallRatio*tRatio)+(sRecallRatio*sRatio);
            mnCnt=min([tSize sSize]);
            finalTruePositive=ratio*mnCnt;
            f=Clusters.F_measure(finalTruePositive, tSize, sSize);
            p=finalTruePositive/sSize;
            r=finalTruePositive/tSize;
        end

        function c=Concordance(teachBinCnts, studBinCnts)
            tSize=sum(teachBinCnts);
            sSize=sum(studBinCnts);
            if tSize == 0 || sSize == 0
                c = 0;
                return;
            end
            teachBins=studBinCnts>teachBinCnts&teachBinCnts>0;
            studBins=studBinCnts<=teachBinCnts&studBinCnts>0;
            truePositive=sum(teachBinCnts(teachBins))+sum(studBinCnts(studBins));
            
            c = truePositive/(tSize + sSize - truePositive);
        end
        
        function [means, teachPtrs, studPtrs, teachWeights, studWeights, ...
                minBinSize]=CreateOld(teachData, studData, minBinSize, ...
                doFixedSplits, cacheFile)
            if nargin < 4
                doFixedSplits=false;
            end
            doNotDuplicate=isequal(studData, teachData);
            if doNotDuplicate
                data=teachData;
                studData=[];
                [N, m] = size(data);
                minN=N;
            else
                data=[teachData;studData];
                [N, m] = size(data);
                minN=min([size(teachData, 1) size(studData,1)]);
            end
            currentBins=cell(1,1);
            firstEntry=cell(1);
            firstEntry{1}=data; %first set is 1 bin-->the original data
            currentBins{1}=firstEntry;
            currentBinPtrs{1}={(1:N)'};
            if nargin<3 || isempty(minBinSize)
                minBinSize=floor(2*log(minN));
            else
                minBinSize=floor(minBinSize);
            end

            i=0;
            while true
                i=i+1;
                if doFixedSplits && i >minBinSize
                    break;
                end
                lastBins=currentBins{i};
                lastBinPtrs=currentBinPtrs{i};
                eventsPerBin=size(lastBins{1}, 1);
                if eventsPerBin==1
                    break;
                end
                if eventsPerBin/2<minBinSize
                    if ~doFixedSplits
                        break;
                    end
                end
                nextBins=cell(1,2^i);
                nextBinPtrs=cell(1,2^i);
                binIdx=1;
                numOfBins=length(lastBins);
                for j=1:numOfBins
                    bin=lastBins{j};
                    binPtrs=lastBinPtrs{j};
                    if size(bin,1) == 1
                        nextBins{binIdx} = bin;
                        nextBinPtrs{binIdx}=binPtrs;
                    else
                        variance=var(bin,1);
                        [~,maxCol]=max(variance);
                        [~,ptrs] = sort(bin(:,maxCol));
                        splitPtr = ceil(length(ptrs)/2);
                        nextBins{binIdx}=bin(ptrs(1:splitPtr),:);
                        nextBinPtrs{binIdx}=binPtrs(ptrs(1:splitPtr),:);
                        binIdx=binIdx+1;
                        nextBins{binIdx}=bin(ptrs(splitPtr + 1:length(ptrs)),:);
                        nextBinPtrs{binIdx}=binPtrs(ptrs(splitPtr+1:length(ptrs)),:);
                    end
                    binIdx=binIdx+1;
                end
                currentBins{i + 1} = nextBins(1:binIdx-1);
                currentBinPtrs{i + 1} = nextBinPtrs(1:binIdx-1); 
            end
            N1 = size(teachData,1);
            teachPtrs=zeros(1, N1);
            finalBins = currentBins{end};
            finalBinPtrs = currentBinPtrs{end};
            numBins=length(finalBins);
            means=zeros(numBins,m);
            if nargout>3
                N2 = size(studData,1);
                studPtrs=zeros(1, N2);
                teachWeights = zeros(numBins,1);
                studWeights = zeros(numBins,1);
                for i=1:numBins
                    bins=finalBins{i};
                    means(i,:)=mean(bins,1);
                    binPtrs=finalBinPtrs{i};
                    p=binPtrs<=N1;
                    teachPtrs(binPtrs(p))=i;
                    studPtrs(binPtrs(~p)-N1)=i;
                    teachWeights(i)=sum(p)/N1;
                    studWeights(i)=sum(~p)/N2;
                end
            else
                if doNotDuplicate
                    for i=1:numBins
                        bins=finalBins{i};
                        means(i,:)=mean(bins,1);
                        binPtrs=finalBinPtrs{i};
                        p=binPtrs<=N1;
                        teachPtrs(binPtrs(p))=i;
                    end
                else
                    N2 = size(studData,1);
                    studPtrs=zeros(1, N2);
                    for i=1:numBins
                        bins=finalBins{i};
                        means(i,:)=mean(bins,1);
                        binPtrs=finalBinPtrs{i};
                        p=binPtrs<=N1;
                        teachPtrs(binPtrs(p))=i;
                        studPtrs(binPtrs(~p)-N1)=i;
                    end
                end
            end
           if AdaptiveBins.DEBUG_LEVEL>0
               teachCnt=floor(mean(teachWeights(teachWeights>0))*N1);
               studCnt=floor(mean(studWeights(studWeights>0))*N2);
               fprintf('this=%d/%d %s%% and stud=%d/%d %s%%\n', ...
                   teachCnt, N1, num2str(round(teachCnt/N1*100, 3)),...
                   studCnt, N2, num2str(round(studCnt/N2*100, 3)));
           end
           if doNotDuplicate
               studPtrs=teachPtrs;
               if nargout>3
                   studWeights=teachWeights;
               end
           end
           if nargin>4 && ~isempty(cacheFile)
               binSize=minBinSize;
               save(cacheFile, 'means', 'teachPtrs', 'studPtrs', 'binSize'); 
           end
        end

        function [means, teachPtrs, studPtrs, teachWeights, studWeights, ...
                minBinSize]=Create(teachData, studData, minBinSize, ...
                doFixedSplits, cacheFile)
            if nargin < 4
                doFixedSplits=false;
            end

            usingTwoDataSets=(nargin>1)&&~isempty(studData)&&~isequal(studData, teachData);
            if usingTwoDataSets
                data=[teachData; studData]; 
                sizeOfSmallerSet=min([size(teachData, 1), size(studData,1)]);
            else
                data=teachData;
                sizeOfSmallerSet=size(teachData, 1);
            end
            [nDataPoints, nParameters] = size(data);

            %It saves time to put both data and pointers in the same
            %array, since the array manipulations throughout the algorithm
            %are identical.
            currentBinsAndPtrs=[data (1:nDataPoints)'];
            eventsPerBin=nDataPoints;
            numOfBins=1;
            if nargin<3 || isempty(minBinSize)
                minBinSize=floor(2*log(sizeOfSmallerSet));
            else
                minBinSize=floor(minBinSize);
            end

            splitNumber=0;
            while true
                splitNumber=splitNumber+1;

                if isTimeToStopSplitting()
                    break;
                end

                rowsForSplittingBins=sortBinsByHighestVarianceParameter();
                currentBinsAndPtrs=splitBinsAccordingToRows();

                eventsPerBin=size(currentBinsAndPtrs, 1);
                numOfBins=size(currentBinsAndPtrs,3);
            end

            means=reshape(mean(currentBinsAndPtrs,'omitnan'),nParameters+1,numOfBins)';
            means(:,end)=[];
            finalPtrs=reshape(currentBinsAndPtrs(:,end,:),eventsPerBin,numOfBins);
            N1 = size(teachData,1);
            teachPtrs=zeros(1,N1);
            teachWeights=zeros(numOfBins,1);

            if usingTwoDataSets
                N2 = size(studData,1);
                studPtrs=zeros(1,N2);
                if nargout > 3
                    studWeights=zeros(numOfBins,1);
                end
                for bin=1:numOfBins
                    ptrsInBin=finalPtrs(:,bin);
                    ptrsInBin(isnan(ptrsInBin))=[];
                    isTeacherData=ptrsInBin<=N1;
                    teachPtrs(ptrsInBin(isTeacherData))=bin;
                    studPtrs(ptrsInBin(~isTeacherData)-N1)=bin;
                    if nargout > 3
                        teachWeights(bin)=sum(isTeacherData)/N1;
                        studWeights(bin)=sum(~isTeacherData)/N2;
                    end
                end
            else
                for bin=1:numOfBins
                    ptrsInBin=finalPtrs(:,bin);
                    ptrsInBin(isnan(ptrsInBin))=[];
                    teachPtrs(ptrsInBin)=bin;
                    teachWeights(bin)=length(ptrsInBin)/N1;
                end
                studPtrs=teachPtrs;
                studWeights=teachWeights;
            end
                
           if AdaptiveBins.DEBUG_LEVEL>0
               teachCnt=floor(mean(teachWeights(teachWeights>0))*N1);
               studCnt=floor(mean(studWeights(studWeights>0))*N2);
               fprintf('this=%d/%d %s%% and stud=%d/%d %s%%\n', ...
                   teachCnt, N1, num2str(round(teachCnt/N1*100, 3)),...
                   studCnt, N2, num2str(round(studCnt/N2*100, 3)));
           end
           if nargin>4 && ~isempty(cacheFile)
               binSize=minBinSize;
               save(cacheFile, 'means', 'teachPtrs', 'studPtrs', 'binSize'); 
           end

            function stopNow=isTimeToStopSplitting()
                stopNow = (doFixedSplits && splitNumber > minBinSize) ...
                    || (~doFixedSplits && eventsPerBin/2<minBinSize) ...
                    || eventsPerBin==1;
            end

            function rows=sortBinsByHighestVarianceParameter()
                                variances=var(currentBinsAndPtrs,1,'omitnan');
                variances=reshape(variances,nParameters+1,numOfBins);
                [~,maxCols]=max(variances(1:end-1,:));

                arrayForSorting=zeros(eventsPerBin,numOfBins);
                for bin_ = 1:numOfBins
                    arrayForSorting(:,bin_)=currentBinsAndPtrs(:,maxCols(bin_),bin_);
                end
                [~,rows]=sort(arrayForSorting);
            end

            function newBinsAndPtrs=splitBinsAccordingToRows()
                newEventsPerBin=ceil(eventsPerBin/2);
                newBinsAndPtrs=zeros(newEventsPerBin,nParameters+1,2*numOfBins);
                if mod(eventsPerBin,2)==0
                    for bin_ = 1:numOfBins
                        newBinsAndPtrs(:,:,2*bin_-1)=currentBinsAndPtrs(rowsForSplittingBins(1:newEventsPerBin,bin_),:,bin_);
                        newBinsAndPtrs(:,:,2*bin_)=currentBinsAndPtrs(rowsForSplittingBins((newEventsPerBin+1):end,bin_),:,bin_);
                    end
                else
                    for bin_ = 1:numOfBins
                        newBinsAndPtrs(:,:,2*bin_-1)=currentBinsAndPtrs(rowsForSplittingBins(1:newEventsPerBin,bin_),:,bin_);
                        newBinsAndPtrs(:,:,2*bin_)=currentBinsAndPtrs(rowsForSplittingBins(newEventsPerBin:end,bin_),:,bin_);
                    end
                    newBinsAndPtrs(end,:,1:2:end)=NaN;
                end
            end
        end
       
        function lvl=DEBUG_LEVEL
            lvl=0;
        end
        
        function emd=Emd(tData, sData, bins)
            try
                tN=size(tData, 1);
                sN=size(sData, 1);
                mx=max([tN sN]);
                minEvents=floor(2*log(mx));
                numberOfBins=floor(mx/minEvents);
                if numberOfBins<32
                    bins=5;
                elseif numberOfBins<64
                    bins=6;
                elseif numberOfBins<128
                    bins=7;
                elseif nargin<3
                    bins=8;
                end
                [tM, ~, ~, tW]=AdaptiveBins.Create(tData, tData, bins, true);
                [sM, ~, ~, sW]=AdaptiveBins.Create(sData, sData, bins, true);
                [~, emd]=emd_flow(tM, sM, tW, sW, @gdf);
            catch ex
                ex.getReport
                emd=nan;
            end
        end
        
         function bins=MeansWeightsPtrs(data)
             bins=SuhProbabilityBins(data);
         end
        
         function bins=Fit(data, binObj, ptrs)
             space='CityBlock';
             [~, bins]=pdist2(binObj.means, data, space, ...
                 'Smallest', 1);
             if nargin>2
                 difs=abs(bins-ptrs);
                 lgcl=difs>0;                 
                 R=size(binObj.means,1);
                 if R>100
                     mn=500;
                 else
                     mn=R/2;
                 end
                 [~, bins4]=pdist2(binObj.means, data(lgcl,:), ...
                     space, 'Smallest', mn);%top 100
                 N=sum(lgcl);
                 idxs=find(lgcl);
                 trueBinDists=zeros(N, 2);
                 for i=1:N
                     idx=idxs(i);
                     b=bins4(:, i);
                     ptr=ptrs(idx);
                     D=pdist2(data(idx,:), binObj.means(ptr,:), space);
                     idx2=find(b==ptr, 1);
                     if ~isempty(idx2)
                         trueBinDists(i,1)=idx2;
                     end
                     trueBinDists(i,2)=D;
                 end
             end
         end
    end 
end
