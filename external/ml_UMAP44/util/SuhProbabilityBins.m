classdef SuhProbabilityBins < handle
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(SetAccess=private)
        means;
        ptrs;
        weights;
        originalData;
        originalLabels;
        lastDecompress;
        retainData=false;
    end
    
    methods
        function this=SuhProbabilityBins(data, retainOriginalData)
            MIN_BINS=8192;
            MIN_EVENTS_PER_BIN=5;
            MAX_EVENTS_PER_BIN=34;
            N=size(data, 1);
            eventsPerBin=floor(2*log(N));
            numberOfBins=floor(N/eventsPerBin);
            if numberOfBins<MIN_BINS
                eventsPerBin=floor(N/MIN_BINS);
            end
            if eventsPerBin<MIN_EVENTS_PER_BIN
                eventsPerBin=MIN_EVENTS_PER_BIN;
            elseif eventsPerBin>MAX_EVENTS_PER_BIN
                eventsPerBin=MAX_EVENTS_PER_BIN;
            end
            if numberOfBins>2^14 %16384
                eventsPerBin=MAX_EVENTS_PER_BIN;
            end
            if nargin>1  && retainOriginalData
                this.originalData=data;
                this.retainData=true;
            end
            [this.means, this.ptrs, ~, this.weights]=...
                AdaptiveBins.Create(data, data, eventsPerBin, false);
        end
        
        function out=compress(this)
            out=this.means;
        end
        
        function out=decompress(this, data, jiggle)
            if nargin<3
                jiggle=[];
            end
            if size(data, 1) ~= size(this.means, 1)
                data=data';
            end
            out=data(this.ptrs,:);
            if this.retainData
                this.lastDecompress=data;
            end
            if ~isempty(jiggle)
                if nargin<4
                    jiggle=.005;
                end
                [R,C]=size(out);
                nudge=((jiggle/2)-...
                    (rand(R,C).*jiggle)).*range(out); 
               out=out+nudge;
            end
        end
        
        function out=fit(this, labels)
            N=size(this.means, 1);
            assert(max(this.ptrs)<=N, 'Expect bin pointer can be higher than bin');
            if size(this.ptrs,2)==1
                assert(size(this.ptrs,1)>N, 'Expect more bin pointers than bins');
            elseif size(this.ptrs,1)==1
                assert(size(this.ptrs,2)>N, 'Expect more bin pointers than bins');
            else
                assert(false, 'Expect 1 column for bin pointers');
            end
            if length(this.ptrs) ~= size(labels,1)
                labels=labels';
            end
            C=size(labels,2);
            out=zeros(N,C);
            out(this.ptrs,:)=labels;
            if this.retainData
                this.originalLabels=labels;
            end
        end
        
    end
    
    methods(Static)
        function eventsPerBin=EventsPerBin(N)
            MIN_BINS=8192;
            MIN_EVENTS_PER_BIN=5;
            MAX_EVENTS_PER_BIN=34;
            eventsPerBin=floor(2*log(N));
            numberOfBins=floor(N/eventsPerBin);
            if numberOfBins<MIN_BINS
                eventsPerBin=floor(N/MIN_BINS);
            end
            if eventsPerBin<MIN_EVENTS_PER_BIN
                eventsPerBin=MIN_EVENTS_PER_BIN;
            elseif eventsPerBin>MAX_EVENTS_PER_BIN
                eventsPerBin=MAX_EVENTS_PER_BIN;
            end
            if numberOfBins>2^14 %16384
                eventsPerBin=MAX_EVENTS_PER_BIN;
            end
        end
        
        function ab=Bins(tData, sData)
            if nargin<2
                sData=tData;
            end
            if isequal(tData, sData)
                N=size(tData, 1);
            else
                N=min(size(tData,1), size(sData,1));
            end
            eventsPerBin=SuhProbabilityBins.EventsPerBin(N);
            ab=AdaptiveBins(tData, sData, eventsPerBin, false);
        end
    end
end