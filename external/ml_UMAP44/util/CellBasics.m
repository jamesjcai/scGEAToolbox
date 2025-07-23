classdef CellBasics
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    methods(Static)
        function tabLines=ToTabLines(data, labels, removeHtml)
            if nargin<3
                removeHtml=true;
            end
            if removeHtml
                J=edu.stanford.facs.swing.Basics;
            end
            tabLines=java.util.ArrayList;
            [rows,cols]=size(data);
            if nargin>1
                data=[labels; data];
            end
            for row=1:rows
                line='';
                for col=1:cols
                    value=data{row, col};
                    if removeHtml && ischar(value)
                        value=J.RemoveXml(value);
                    end
                    if isnumeric(value)
                        value=num2str(value);
                    end
                    line=sprintf('%s%s\t', line, value);                    
                end
                tabLines.add(line);
            end
        end
        
        function idx=IndexOf(items, item)
            idx=0;
            N=length(items);
            for i=1:N
                if isequal(items{i}, item)
                    idx=i;
                    break;
                end
            end
        end
        
        function groups=GetRowGroups(sortedRows)
            groups={};
            R=size(sortedRows, 1);
            if R>0
                lastStartingR=1;
                for r=2:R
                    if sortedRows(r,1)>sortedRows(r-1,1)+1
                        groups{end+1}=sortedRows(lastStartingR:r-1,:);
                        lastStartingR=r;
                    end
                end
            end
            groups{end+1}=sortedRows(lastStartingR:r,:);
        end
        
        function rowIdxs=Find(c, colIdxs, fnc, first, rowStartIdx)
            if nargin<5
                rowStartIdx=1;
                if nargin<4
                    first=true;
                end
            end
            N=size(c,1);
            C=size(c,2);
            colIdxsN=length(colIdxs);
            rowIdxs=[];
            for i=rowStartIdx:N
                for j=1:colIdxsN
                    colIdx=colIdxs(j);
                    if colIdx>=1 && colIdx<=C
                        o=c{i,colIdx};
                        if fnc(o, i, colIdx)
                            if colIdxsN==1
                                rowIdxs(end+1)=i;
                            else
                                rowIdxs(end+1,:)=[i colIdx];
                            end
                            if first
                                return;
                            end
                        end
                    end
                end
            end
        end
        
        function yes=equals(o, arg)
            yes=o==arg;
        end
        
        function [found, idx]=Contains(cell, object)
            found=false;
            N_=length(cell);
            for idx=1:N_
                if isequal(object, cell{idx})
                    found=true;
                    return;
                end
            end       
            idx=0;
        end

        function c=Java(collectionOrArray)
            try
                N=collectionOrArray.size;
                c=cell(1,N);
                it=collectionOrArray.iterator;
                for i=1:N
                    c{i}=it.next;
                end
            catch 
                N=length(collectionOrArray);
                c=cell(1,N);
                for i=1:N
                    c{i}=collectionOrArray(i);
                end
                
            end
        end

        function u=UniqueNumbers(cellNums)
            if isnumeric(cellNums)
                u=unique(cellNums);
            else
                u=[];
                N=length(cellNums);
                for i=1:N
                    u=unique([u cellNums{i}']);
                end
            end
        end
        
        function [in, I]=Order(canon, in)
            N1=length(canon);
            N2=length(in);
            if N1==0
                I=1:N2;
                return;
            end
            foundIdxs=cell(1,N1);
            newIdxs=[];
            for i=1:N2
                item=in{i};
                found=false;
                for j=1:N1
                    if isequal(item, canon{j})
                        if isempty(foundIdxs{j})
                            foundIdxs{j}={i};
                        else
                            foundIdxs{j}=[foundIdxs{j} i];
                        end
                        found=true;
                        break;
                    end
                end
                if ~found
                    newIdxs(end+1)=i;
                end
            end
            I=zeros(1,N2);
            N3=length(foundIdxs);
            idx=1;
            for i=1:N3
                N4=length(foundIdxs{i});
                for j=1:N4
                    I(idx)=foundIdxs{i}{j};
                    idx=idx+1;
                end
            end
            N5=length(newIdxs);
            for i=1:N5
                I(idx)=newIdxs(i);
                idx=idx+1;
            end
            in=in(I);
        end
    end
end
