classdef LabelBasics < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    methods(Static)
        function [labels, loss]=RemoveOverlap(labels)
            nCols=size(labels,2);
            if nCols>1
                y=labels(:,1);
                nCols=size(labels,2);
                
                for col=2:nCols
                    ll=y==0;
                    more=labels(ll, col);
                    y(ll)=more;
                end
                if nargout>1
                    lbls=unique(labels);
                    lbls=lbls(lbls>0);
                    h1=MatBasics.HistCounts(sort(labels(:)), lbls);
                    h2=MatBasics.HistCounts(sort(y), lbls);
                    loss=h2-h1;
                end
                labels=y;
            end
        end
        
        % assumes AutoGate modules are on path
        function map=GetLabelMap(gt, labels)
            ids=unique(labels);
            [names, N, strIds]=GatingTree.GetUniqueNames(gt.tp, ids);
            map=java.util.Properties;
            hl=gt.highlighter;
            for i=1:N
                dfltClr=Gui.HslColor(i,N);
                if ids(i)~=0
                    id=strIds{i};
                    name=names{i};
                    if isempty(name)
                        name=['label=' id];
                    end
                    if gt.tp.hasChildren(id)
                        name=[name ' (non-specific)'];
                    end
                    map.put(java.lang.String(id), name);
                    clr=num2str(floor(hl.getColor(id, dfltClr)*256));
                    map.put([id '.color'], clr);
                end
            end
        end
        
    
        function [key, keyColor, keyTraining, keyTest]=Keys(lbl)
            key=num2str(lbl); % prevent java.lang.Character
            keyColor=[key '.color'];
            keyTraining=[key '.trainingFrequency'];
            keyTest=[key '.testFrequency'];
            key=java.lang.String(key);
        end
        
        function Frequency(lbls, lblMap, training, match, how)
            if nargin<5
                how=-1; %contains, 0=eq, 1 = startsWith
                if nargin<4
                    match='';
                    if nargin<3
                        training=true;
                    end
                end
            end
            [trainFreq, ~]=LabelBasics.HasFrequencies(lblMap);
            total=length(lbls);
            tab=sprintf('\t');
            u=unique(lbls);
            nLbls=length(u);
            frequencies=MatBasics.HistCounts(sort(lbls), u);
            [sortFrequencies,II]=sort(frequencies, 'descend');
            if nargin>3
                if ~isempty(training) && training
                    disp([num2str(nLbls) ' training set labels']);
                else
                    disp([num2str(nLbls) ' test set labels']);
                end
            end
            for i=1:nLbls
                lbl=u(II(i));
                freq=sortFrequencies(i);
                pFreq=String.encodePercent(freq, total, 1);
                [key, ~, keyTraining, keyTest]=LabelBasics.Keys(lbl);
                if ~isempty(training)
                    if training
                        lblMap.put(keyTraining, pFreq);
                    else
                        lblMap.put(keyTest, pFreq);
                    end
                end
                if nargin>3
                    name=lblMap.get(key);
                    if trainFreq && ~training
                        start=['From ' lblMap.get(keyTraining) ' to '];
                    else
                        start='';
                    end
                    if ~isempty(match)
                        start=[tab start];
                        if how==-1
                            if String.Contains(name, match)
                                start=['* ' start];
                            end
                        elseif how==1
                            if String.StartsWith(name, match)
                                start=['* ' start];
                            end
                        else
                            if strcmp(name, match)
                                start=['* ' start];
                            end
                        end
                    end
                    sFreq=String.encodeInteger(freq);
                    fprintf('%s%s=%s events for #%d="%s"\n', start, ...
                        pFreq, sFreq, lbl, name);
                end
            end
            if nargin>3
                if ~isempty(training) && training
                    disp([num2str(nLbls) ' training set labels']);
                else
                    disp([num2str(nLbls) ' test set labels']);
                end
            end
            if ~isempty(training)
                if training
                    lblMap.put('hasTrainingFrequencies', 'yes');
                else
                    lblMap.put('hasTestFrequencies', 'yes');
                end
            end
        end
        
        function [training, test]=HasFrequencies(lblMap)
            if isempty(lblMap)
                training=false;
                test=false;
            else
                training=strcmpi('yes', lblMap.get('hasTrainingFrequencies'));
                test=strcmpi('yes', lblMap.get('hasTestFrequencies'));
            end
        end
        
        function [addTrainingHtml, sup1, sup2, trStart, trEnd]=...
                AddTrainingHtml(lblMap, needHtml)
            addTrainingHtml=false;
            if needHtml
                sup1=BasicMap.Global.supStart;
                sup2=BasicMap.Global.supEnd;
                if ~isempty(lblMap)
                    trStart=[sup1 '<font color="#11DDAA"><b> training '];
                    trEnd=['</b></font>' sup2];
                    [trainFreq, testFreq]=LabelBasics.HasFrequencies(lblMap);
                    addTrainingHtml=trainFreq&&testFreq;
                else
                    trStart=[]; trEnd=[];
                end
            else
                sup1=[]; sup2=[]; trStart=[]; trEnd=[];
            end

        end
        
        function[data, columnNames1]=Merge2Samples(...
                fileName, sample1, label1, sample2, label2, label3)
            [data1, columnNames1]=File.ReadCsv(sample1);
            [data2, columnNames2]=File.ReadCsv(sample2);
            if ~isequal(columnNames1(1:end-1), columnNames2(1:end-1))
                msgWarning('Training & test set labels do not match');
                return;
            end
            lblMap2=loadLabels(label2);
            lblMap1=loadLabels(label1);
            if isequal(label3, label1)
                lbls=data2(:,end);
                reLbls=LabelBasics.RelabelIfNeedBe(lbls, lblMap1, lblMap2);
                if size(reLbls,2)~=size(lbls,2)
                    data=[];
                    return;
                end
                data2(:,end)=reLbls;
            else
                lblMap3=loadLabels(label3);
                lbls=data1(:,end);
                reLbls=LabelBasics.RelabelIfNeedBe(lbls, lblMap3, lblMap1);
                if size(reLbls,2)~=size(lbls,2)
                    data=[];
                    return;
                end
                data1(:,end)=reLbls;
                lbls=data2(:,end);
                reLbls=LabelBasics.RelabelIfNeedBe(lbls, lblMap3, lblMap2);
                if size(reLbls,2)~=size(lbls,2)
                    data=[];
                    return;
                end
                data2(:,end)=reLbls;
            end
            data=[data1;data2];
            if ~isempty(fileName)
                try
                    fu=edu.stanford.facs.wizard.FcsUtil(fileName);
                catch ex
                    msgError('Cannot load java edu.stanford.facs.wizard.FcsUtil' );
                end
                problem=fu.createTextFile(fileName, [], data, [], columnNames1);
                copyfile(label1, File.SwitchExtension(fileName, '.properties')) 
                if ~isempty(problem)
                    msgError(problem, 12);
                end
            end
            
            function lblMap=loadLabels(lblFile)
                try
                    lblMap=java.util.Properties;
                    lblMap.load(java.io.FileInputStream(lblFile));
                catch ex
                    ex.getReport
                    lblMap=[];
                end
            end
        end
        
        function lbls=RelabelIfNeedBe(lbls, trainingMap, testMap)
            if isempty(trainingMap) || trainingMap.size()==0
                msgWarning('Training set map is empty');
                lbls=[];
                return;
            end
            if isempty(testMap) || testMap.size()==0
                msgWarning(Html.WrapHr(['No classification labels '...
                    '<br>for test set.']));
                lbls=[];
                return;
            end
            u=unique(lbls);
            N=length(u);
            trainingIdByName=LabelBasics.IdByName(trainingMap);
            for i=1:N
                lbl=u(i);
                if lbl ~= 0
                    key=java.lang.String(num2str(lbl));
                    if ~trainingMap.containsKey(key)
                        name=testMap.get(key);
                        if isempty(name)
                           warning(['Test set properties lack'...
                               ' name for label "' char(key) '"']); 
                        else
                            newLbl=trainingIdByName.get(name);
                            if ~isempty(newLbl)
                                newLbl=str2double(newLbl);
                                lbls(lbls==lbl)=newLbl;
                            else
                                warning(['Training set properties lack'...
                                    ' name for label "' char(key) '"'...
                                    'named "' name '"']);
                            end
                        end
                    end
                end
            end            
        end
        
        function map=IdByName(inMap)
            map=java.util.TreeMap;
            it=inMap.keySet.iterator;
            while it.hasNext
                key=char(it.next);
                if ~endsWith(key, '.color')
                    name=inMap.get(key);
                    if map.containsKey(key)
                        warning(['Duplicate use of ' name]);
                    else
                        map.put(java.lang.String(name), key);
                    end
                end
            end
        end
        
         
        function counts = DiscreteCount(x, labels)
            if isempty(labels)
                counts = [];
                return;
            end
            if ~any(isinf(labels))
                labels(end+1) = inf;
            end
            counts = histcounts(x, labels);
        end
        
        function [names, clrs, lbls]=GetNamesColorsInLabelOrder(...
                lbls, lblMap, minFrequency)
            lbls(lbls<0)=0;
            ids_=unique(lbls);
            N_=length(ids_);
            cnts_ = LabelBasics.DiscreteCount(lbls, ids_);
            if nargin>2
                isSigId = cnts_>=minFrequency;
            else
                isSigId=true(1, N_);
            end
            if all(size(ids_)== size(isSigId))
                ids_=ids_';
            end
            N_Names = sum(isSigId & (ids_ > 0)');
            names=cell(1,N_Names);
            clrs=zeros(N_Names,3);
            sig_idx=0;
            for i=1:N_
                id=ids_(i);
                if id>0
                    key=num2str(id);
                    name=lblMap.get(java.lang.String(key));
                    if isempty(name)
                        name=['Subset #' key];
                    end
                    if isSigId(i)
                        sig_idx=sig_idx+1;
                        key=num2str(id);
                        names{sig_idx}=name;
                        try
                            clrs(sig_idx,:)=str2num(...
                                lblMap.get([key '.color']))/256; %#ok<ST2NM>
                        catch 
                            clrs(sig_idx,:)=[.1 .1 .1];
                        end
                    else
                        if ~strcmpi(this.verbose, 'none')
                            warning('Ignoring "%s" since it has %d events...', ...
                                name, cnts_(i));
                        end
                        lbls(lbls==id)=0;
                    end
                else
                    lbls(lbls==id)=0;
                end
            end
        end
       
        function [map, halt, args]=GetOrBuildLblMap(lbls, args)
            halt=false;
            map=[];
            if isempty(args.label_file)
                warning(['label_column without label_file '...
                    'to match/supervise, will use default names/colors']);
                args.buildLabelMap=true;
            end
            if args.buildLabelMap
            else
                if ~exist(args.label_file, 'file')
                    label_file=WebDownload.GetExampleIfMissing(args.label_file);
                    if exist(label_file, 'file')
                        args.label_file=label_file;
                    end
                end
                if exist(args.label_file, 'file')
                    map=File.ReadProperties( args.label_file);
                    if isempty(map)
                        problem='load';
                    end
                elseif ~isempty(args.label_file)
                    problem='find';
                end
                if isempty(map)
                    globals=BasicMap.Global;
                    if askYesOrNo(['<html>Cannot ' problem ' the '...
                            ' label file <br><br>"<b>' globals.smallStart ...
                            args.label_file globals.smallEnd '</b>"<br><br>'...
                            '<center>Use default names & black/white colors?</center>'...
                            '<hr></html>'], 'Error', 'north west', true)
                        args.buildLabelMap=true;
                    else
                        halt=true;
                    end
                end
            end
            if args.buildLabelMap
                map=java.util.Properties;
                u=unique(lbls)';
                nU=length(u);
                if nU/args.n_rows > .2
                    if ~acceptTooManyLabels(nU, args.n_rows)
                        halt=true;
                        return;
                    end
                end
                for i=1:nU
                    key=num2str(u(i));
                    map.put(java.lang.String(key), ['Subset #' key]);
                    map.put([key '.color'], num2str(Gui.HslColor(i, nU)));
                end
            end
            if args.color_defaults
                ColorsByName.Override(map, args.color_file, beQuiet);
            end
            
            function ok=acceptTooManyLabels(nU, nRows)
                ok=true;
                txt=sprintf(['You have %s unique labels...<br>'...
                    'This is %s of the actual data rows...'], ...
                    String.encodeInteger(nU), String.encodePercent(...
                    nU/nRows, 1, 1));
                if ~askYesOrNo(Html.WrapHr(...
                        sprintf(['Interesting ....%s'...
                        '<br><br>So this then will be very SLOW ...'...
                        '<br><br><b>Continue</b>????'], ...
                        txt)))
                    ok=false;
                    return;
                end
            end
        end

        function [ok, cancelled]=Confirm(labels, limit, askToTreatAsData)
            if nargin<3
                askToTreatAsData=true;
            end
            cancelled=false;
            u=unique(labels);
            warningTxt='';
            nFound=length(u);
            strN=String.encodeInteger(nFound);
            total=length(labels);
            if limit<1
                if nFound/total>limit
                    warningTxt=sprintf(['<u>%s</u><br>...yet %s ' ...
                        '<i>is <u><b>%s</b></u> of %s</i><hr>'], ...
                        String.encodePercent(limit), strN, ...
                        String.encodePercent(nFound/total),  ...
                        String.encodeInteger(total));
                end
            else
                if nFound>limit
                    warningTxt=sprintf('<u>%s</u>!<hr>',...
                        String.encodeInteger(limit));
                end                
            end
            if ~isempty(warningTxt)
                question=[strN ' classification labels have been'...
                    ' <br>found in ' String.encodeInteger(total) ...
                    ' matrix rows!<br>...BUT the limit is ' warningTxt];
                if askToTreatAsData
                    [a, cancelled]=Gui.Ask(Html.WrapHr(question), {...
                        'These are labels', ...
                        'Treat as data',...
                        'STOP'}, 'LabelBasics.Confirm');
                    if ~cancelled
                        ok=a==1;
                        cancelled=a==3;
                    else
                        ok=false;
                    end

                else
                    ok=askYesOrNo(Html.WrapC([question '<br><br><b>'...
                        '<font color="red">Continue?</font></b>']));
                end
            else
                ok=true;
            end
        end

        function table=CompressTable(table, ...
                dataSetFactor, classFactor, columnNames)
            inData=table{:,1:end-1};
            labels=table{:,end};
            [inData, labels]=LabelBasics.Compress(...
                inData, labels, dataSetFactor, ...
                classFactor);
            table=array2table([inData labels], ...
                'VariableNames', columnNames);
        end

        function [data, labels]=Compress(data, labels, ...
                dataSetFactor, classFactor)
            if nargin<4
                classFactor=0;
            end
            R=size(data, 1);
            if dataSetFactor==1
                return;
            end
            if floor(dataSetFactor)~=dataSetFactor
                sz=floor(dataSetFactor*R);
            else
                sz=dataSetFactor;
            end
            if sz>=R
                warning('%s factor produces %d and data set has %d rows',...
                    String.encodeRounded(dataSetFactor,2), sz, sz);
                return;
            end
            r=randperm(R);
            r=r(1:sz);
            if ~isempty(labels)
                if classFactor>0
                    uOriginal=unique(labels)';
                    cntsOriginal=MatBasics.HistCounts(labels, uOriginal);
                    smallestSubset=min(cntsOriginal);
                    isRatio=floor(classFactor)~=classFactor;
                    if isRatio
                        minimum=classFactor*smallestSubset;
                    else
                        minimum=classFactor;
                        if any(minimum>cntsOriginal)
                            warning(...
                                ['BEFORE compression %d is > '...
                                'than %d subset(s):'...
                                '  labels=[%s]. Counts=[%s]'], ...
                                minimum, sum(minimum>cntsOriginal), ...
                                MatBasics.toString(...
                                uOriginal(minimum>cntsOriginal)), ...
                                MatBasics.toString(...
                                cntsOriginal(minimum>cntsOriginal))...
                                );
                        end
                    end
                    temp=labels(r);
                    uCompressed=unique(temp)';
                    cntsCompressed=MatBasics.HistCounts(temp, uCompressed);
                    missing=uOriginal(~ismember(uOriginal, uCompressed));
                    tooSmallR=uCompressed(minimum>cntsCompressed);
                    needToIncrease=[missing tooSmallR];
                    nTooSmallLabels=length(needToIncrease);
                    if nTooSmallLabels>0
                        r=r(~ismember(temp, needToIncrease));
                        for i=1:nTooSmallLabels
                            label=needToIncrease(i);
                            idxs=find(labels==label)';
                            nIdxs=length(idxs);
                            if nIdxs<=minimum
                                r=[r idxs];
                            else
                                r2=randperm(nIdxs);
                                r=[r idxs(r2(1:minimum))];
                            end
                        end
                    end
                    data=data(r,:);
                    labels=labels(r);
                    if nTooSmallLabels>0
                        cntsCompressedScaledUp=MatBasics.HistCounts(labels, uCompressed);
                        originalTooSmall=uOriginal(minimum>cntsOriginal);
                        newTooSmall=~ismember(needToIncrease, originalTooSmall);
                        if any(newTooSmall)
                            % make nice report for warning
                            newTooSmallLabels=needToIncrease(newTooSmall);
                            N2=length(newTooSmallLabels);
                            s='';
                            for j=1:N2
                                label=newTooSmallLabels(j);
                                cnt1=cntsOriginal(find(uOriginal==label,1));
                                idx2=find(uCompressed==label,1);
                                if idx2==0
                                    cnt2=0;
                                else
                                    cnt2=cntsCompressed(idx2);
                                end
                                s=[s sprintf('[%d %d %d] ', label, cnt1, cnt2)];
                            end
                            warning('labels=[%s]. Counts=[%s]', ...
                                MatBasics.toString(uCompressed), ...
                                MatBasics.toString(cntsCompressedScaledUp));
                            warning(...
                                ['AFTER compression %d subsets '...
                                'were scaled up to %d),'...
                                '  [label before compressed]:  %s'], ...
                                N2, minimum, s); 
                        end
                    end
                else
                    data=data(r,:);
                    labels=labels(r);
                end
            else
                data=data(r,:);
            end
        end
        
        function key=ColorKey(key)
            key=[key '.color'];
        end

        function [mdns, sizes]=Median(data, labels, uniqueLabels, ...
                ignore0, columnNames, labelProps, filterColumns)
            if nargin<4
                ignore0=true;
                if nargin<3
                    uniqueLabels=[];
                    if nargin<2
                        labels=[];
                    end
                end
            end
            if isempty(labels)
                labels=data(:,end);
                data=data(:,1:end-1);
            end
            if isempty(uniqueLabels)
                uniqueLabels=unique(labels);
            end
            C=size(data,2);
            N=length(uniqueLabels);
            mdns=[];
            sizes=[];
            for i=1:N
                u=uniqueLabels(i);
                if ignore0 && u==0
                    continue;
                end
                mdns(end+1,:)=median(data(labels==u, :));
                sizes(end+1)=sum(labels==u);
            end
            if nargin>4 && nargin>5
                if nargin>6
                    N3=length(filterColumns);
                end
                for j=1:C
                    cn=columnNames{j};
                    idx=String.IndexOf(cn, ':');
                    if idx>0
                        columnNames{j}=cn(1:idx-1);
                    end
                end
                iMdn=0;
                for i=1:N
                    u=uniqueLabels(i);
                    if ignore0 && u==0
                        continue;
                    end
                    iMdn=iMdn+1;
                    name=labelProps.get(num2str(u));
                    if isempty(name)
                        name='background';
                    end
                    fprintf('%s ', name);
                    for j=1:C
                        cn=columnNames{j};
                        if nargin>6
                            found=false;
                            for k=1:N3
                                if contains(cn, filterColumns{k})
                                    found=true;
                                    break;
                                end
                            end
                            if ~found
                                continue;
                            end
                        end
                        fprintf('%s=%s ', cn, ...
                            String.encodeRounded(mdns(iMdn,j), 3));
                    end
                    fprintf('\n');
                end
            end
        end

        function [names, clrs, u]=NamesAndColors(lbls, idMap)
            u=unique(lbls)';
            N=length(u);
            names={};
            clrs=[];
            for i=1:N
                id=num2str(u(i));
                name=idMap.get(id);
                if isempty(name)
                    continue;
                end
                names{end+1}=name;
                clr=idMap.get([id '.clr']);
                if isempty(clr)
                    clr=idMap.get([id '.color']);
                    if isempty(clr)
                        clr=[.8 .8 .84];
                    elseif ischar(clr)
                        clr=str2num(clr)/256; %#ok<ST2NM> 
                    end
                end
                clrs(end+1,:)=clr;
            end
        end

        function labelMap=EmptyMap(labels)
            labelMap=JavaProperties;
            u=unique(labels);
            nU=length(u);
            for i=1:nU
                key=num2str(u(i));
                labelMap.set(java.lang.String(key), ['Subset ' key]);
            end
        end
        
        function outMap=InvertMap(inMap)
            outMap=JavaProperties;
            if isa(inMap, 'java.util.Map')
                it=inMap.keySet.iterator;
                while it.hasNext
                    key=char(it.next);
                    value=inMap.get(key);
                    outMap.set(value, key);
                end
            else
                keys=inMap.keys;
                N=length(keys);
                for i=1:N
                    key=keys{i};
                    value=inMap.get(key);
                    outMap.set(value, key);
                end
            end
        end
        
        function numericLabels=ToNumericLabels(stringLabels, map)
            U=unique(stringLabels);
            N=length(U);
            numericLabels=zeros(length(stringLabels), 1);
            map=LabelBasics.InvertMap(map);
            for i=1:N
                try
                    key=U{i};
                    id=str2num(map.get(key));
                    numericLabels(strcmp(key, stringLabels))=id;
                catch
                end
            end
        end
    end


    properties(Constant)
        PROP_EXCLUDE='LabelBasics.Exclusions.';
        PROP_OVERLAP='classOverlap';
        PROP_REMEMBER_OVERLAP=['remember.' LabelBasics.PROP_OVERLAP];
        
    end

    properties(SetAccess=private)
        contextExcludeProperty;
        contextDescription;
        dataPlural;
        dataSingular;
        contextHtml;
        props;
        propsForLabels;
        idPerEvent;
        overFlowStrategy=2;
        doneGates=0;
        nEvents=0;
        descriptions;
        leavesOnly=true; % ignore final classification of an event
                        % if it does NOT occur within a leaf class/gate... 
                        % e.g. ignore B cells not in the B-1, B-2 
                        % or other sub classes
        pidLevels;
        pidSampleRows;
        subsetDescription;
    end
    
    methods
        function setLeafsOnly(this, yes)
            if ~yes
                this.pidLevels=TreeMapOfMany(java.util.TreeMap);
                this.pidSampleRows=java.util.HashMap;
            end
            this.leavesOnly=yes;
        end
        
        function [priorIds, priorIdxs]=getExclusions(this, u)
            priorIdxs=[];
            v=this.props.get(this.contextExcludeProperty);
            priorIds=[];
            if ~isempty(v)
                priorIds=str2num(v); %#ok<ST2NM> 
                nPrior=length(priorIds);
                for i=1:nPrior
                    priorId=priorIds(i);
                    priorIdx=find(u==priorId, 1);
                    if ~isempty(priorIdx)
                        priorIdxs(end+1)=priorIdx;
                    end
                end
            end
        end
        
        function oneColumn=makeOneColumn(this)
            oneColumn=this.idPerEvent(:,1);
            nCols=size(this.idPerEvent,2);
            if nCols>1
                [cnts, idxsBig1st, ids]=MatBasics.HistCounts(this.idPerEvent);
                oneColumn=zeros(size(this.idPerEvent, 1),1);
                N=length(idxsBig1st);
                for i=N:-1:1
                    idx=idxsBig1st(i);
                    id=ids(idx);
                    if id==0
                        continue;
                    end
                    for col=1:nCols
                        l=this.idPerEvent(:,col)==id;
                        if any(l)
                            oneColumn(l&(oneColumn==0))=id;
                            break;
                        end
                    end
                end
                ids2=unique(oneColumn);
                NN=length(ids2);
                if N > NN
                    msgWarning(Html.SprintfHr(['Overlap hides %s !' ...
                        '<br><br>Will try to represent all...'], ...
                        String.Pluralize2('data subset', N-NN)), 12, 'east');
                    ratio=.4;
                    while true
                        %toDo=N-NN;
                        for i=1:N
                            idx=idxsBig1st(i);
                            id=ids(idx);
                            if id==0
                                continue;
                            end
                            cnt=cnts(idx);
                            cnt2=sum(oneColumn==id);
                            mn=ceil(cnt*ratio);
                            if cnt2<mn
                                for col=1:nCols
                                    idxs2=find(this.idPerEvent(:,col)==id);
                                    NN=length(idxs2);
                                    if NN>0
                                        if mn>NN
                                            r=NN;
                                        else
                                            r=mn;
                                        end
                                        idxs2=idxs2(1:r);
                                        oneColumn(idxs2)=id;
                                        break;
                                    end
                                end
                            end
                        end
                        NN=length(unique(oneColumn));
                        nextToDo=N-NN;
                        if nextToDo==0% || nextToDo>=toDo %nothing changed
                            break;
                        end
                        ratio=ratio/2;
                        if ratio<.05
                            break;
                        end
                    end
                    if N>NN
                        msgWarning(Html.SprintfHr(['Sigh ... overlap STILL ' ...
                            'hides %s!'], String.Pluralize2( ...
                            'data subset', N-NN)), 12, 'south east');
                    end
                end
            end
        end
        
        function ids=applyExclusions(this)
            ids=this.idPerEvent;
            v=this.props.get(this.contextExcludeProperty);
            if isempty(v)
                return;
            end
            deSelected=str2num(v); %#ok<ST2NM> 
            N=length(deSelected);
            C=size(this.idPerEvent, 2);
            for i=1:N
                id=deSelected(i);
                for c=1:C
                    l=this.idPerEvent(:,c)==id;
                    if any(l)
                        ids(l, c)=0;
                        break;
                    end
                end
            end
        end
        
        function [classIds, cancelled]=choose(this, make1Column, ...
                parentFig, modal, calledByCheck, askUser)
            if nargin<6
                askUser=true;
                if nargin<5
                    calledByCheck=false;
                    if nargin<4
                        modal=true;
                        if nargin<3
                            parentFig=get(0, 'CurrentFigure');
                            if nargin<2
                                make1Column=1;
                            end
                        end
                    end
                end
            end
            PROP_SUBSET_OVERLAP=[LabelBasics.PROP_OVERLAP ...
                '.' this.contextExcludeProperty];
            classIds=this.applyExclusions;
            
            [~, mapRatios, mapCounts]=...
                MatBasics.GreaterThan(classIds, 1);
            ratio=MatBasics.GetClassifiedOverlap(classIds);
                
            cancelled=false;
            if ratio>0
                oldStandardRatio=this.props.getNumeric(...
                    LabelBasics.PROP_OVERLAP, -1);
                if abs(ratio-oldStandardRatio)<.001
                    if ~calledByCheck && make1Column
                        classIds=this.makeOneColumn;
                    end
                    return;
                end
                oldRatio=this.props.getNumeric(PROP_SUBSET_OVERLAP, -1);
                if oldRatio>=0
                    if abs(ratio-oldRatio)<.001
                        if ~calledByCheck && make1Column
                            classIds=this.makeOneColumn;
                        end
                        return;
                    end
                end
                oldRatio=this.props.getNumeric(...
                    LabelBasics.PROP_OVERLAP, -1);
                if oldRatio>=0
                    if abs(ratio-oldRatio)<.001
                        if ~calledByCheck && make1Column
                            classIds=this.makeOneColumn;
                        end
                        return;
                    end
                end
                strRatioTitle=String.encodePercent(ratio, 1, 1);
                strRatio=String.ToHtml(strRatioTitle);
                u=unique(classIds(:));
                u=u(u>0);
                %sort by GatingTree leaf order
                if ~calledByCheck
                    if askUser
                        check;
                    end
                    if make1Column
                        classIds=this.makeOneColumn;
                    end
                    return;
                end
                U=length(u);
                choices=cell(1,U);
                for i=1:U
                    id=u(i);
                    overlap=mapRatios.get(id);
                    if overlap>0
                        count=mapCounts.get(id);
                        suffix=['<font color=''red''>'...
                            String.encodeInteger(count) ' events ' ...
                            String.encodePercent(overlap, 1, true)...
                            ' overlap</font> ' ...
                            Html.EncodeSort('overlap', count)];
                    else
                        suffix=Html.EncodeSort('overlap', 0);
                    end
                    description=this.descriptions.get(num2str(id));
                    choices{i}=['<html>' description ' ' suffix '</html>'];
                end
                [~, priorIdxs]=this.getExclusions(u);
                while true
                    where='east+';
                    msgObj=struct(...
                        'modal', modal, ...
                        'msg', ['<html><center><b> <font '...
                        'color=''red''>' strRatio '</font> of the '...
                        this.dataPlural ' are classified<br> in <i>2 or more<'...
                        '/i> subsets...</b><br><br>To ensure that '...
                        'each ' this.dataSingular ' occurs in 1 or 0 '...
                        'subsets...<br><b>select</b> subsets '...
                        '<b>to EXCLUDE</b> from the subsets<br>for '...
                        this.contextHtml ...
                        '<hr></center></html>'],...
                        'where', where, 'icon', 'none', 'sortProps', ...
                        BasicMap.Global, 'sortProp', 'EventIds.overlap');
                    if ~isempty(parentFig)
                        msgObj.javaWindow=Gui.JWindow(parentFig);
                    end
                    [idxs, cancelled]=mnuMultiDlg(msgObj,...
                        [strRatioTitle ' subset overlap... '...
                        '...'],  choices, priorIdxs-1, false, true);
                    if length(idxs)==length(choices)
                        msgBox(Html.WrapHr(['You are excluding every '...
                            'subset ??? <br><br>This of course means' ...
                            ' NOTHING will be done ...<br><br>Try again '...
                            '...or select cancel']));
                    else
                        break;
                    end
                end
                if ~cancelled && ~isempty(idxs)                    
                    deSelected=u(idxs)';
                    if size(deSelected, 1)>1
                        v=num2str(deSelected');
                    else
                        v=num2str(deSelected);
                    end
                    this.props.set(this.contextExcludeProperty, v);
                    classIds=this.applyExclusions;
                    check;
                elseif ~cancelled
                    this.props.set(this.contextExcludeProperty, '');
                    check
                end
            end
            if ~calledByCheck
                if make1Column
                    classIds=this.makeOneColumn;
                end
            end

            
            function check
                [classified, classifiedOverlap, nClasses]...
                    =MatBasics.CountClassified(classIds);
                newRatio=classifiedOverlap/classified;
                if newRatio>0
                    subContext=sprintf(['<br>%s events are in ' ...
                        '%s<br>and %s of these events are ' ...
                        '<br>in <i>more than 1</i> %s.<br>'], ...
                        String.encodeK(classified), ...
                        String.Pluralize2(this.subsetDescription, nClasses), ...
                        String.encodeInteger(classifiedOverlap), ...
                        this.subsetDescription);
                    strNewRatio=String.ToHtml(...
                        String.encodePercent(newRatio, 1, 1));
                    noRetry='No, try to improve overlap';
                    yesOk=['<html>Yes <b>' strNewRatio ...
                        '</b> overlap is OK</br></html>'];
                    if newRatio>.24
                        dfltIdx=1;
                    else
                        dfltIdx=2;
                    end
                    cb=Gui.CheckBox('<html>Remember for<br>this context</html>');
                    if strcmp(strRatio, strNewRatio)
                        app=BasicMap.Global;
                        txtMsg=sprintf(['<html><center>For <b>%s</b>'...
                            '%s<br>Is <font color="red">' ...
                            '<b>%s</b></font> overlap <u>acceptable' ...
                            '</u> ???<br>%s(the prior was <b>%s</b>)%s' ...
                            '<hr></center></html>'], ...
                            this.contextHtml, subContext, strNewRatio, ...
                            app.smallStart, strRatio, app.smallEnd);
                    else
                        txtMsg=sprintf(['<html><center>For <b>%s</b>'...
                            '%s<br><br>Is <b>%s</b> overlap' ...
                            ' <u>acceptable</u> ???<hr>'...
                            '</center></html>'], this.contextHtml,...
                            subContext, strNewRatio);
                    end
                    msgObj=struct('msg', txtMsg, 'where', 'center');
                    if ~isempty(parentFig)
                        msgObj.javaWindow=Gui.JWindow(parentFig);
                    end
                    [~,cancelled, answ]=Gui.Ask(msgObj, {noRetry, yesOk}, [], ...
                        'Classification overlap...', dfltIdx, cb);
                    if ~cancelled 
                        if cb.isSelected && isequal(answ, yesOk)
                            this.props.set(PROP_SUBSET_OVERLAP,...
                                num2str(newRatio));
                            if askYesOrNo(['<html>Allow overlap &lt;= '...
                                    strrep(strNewRatio, '&lt;', '')  ...
                                    ' for <b>ALL</b> subsets?<hr>'...
                                    '</html>'], 'New overlap standard?', ...
                                    'center', true, ...
                                    this.PROP_REMEMBER_OVERLAP)
                                this.props.set(this.PROP_OVERLAP, ...
                                    num2str(newRatio));
                            end
                        elseif isequal(answ, noRetry)
                            [classIds, cancelled]=this.choose(...
                                make1Column, parentFig, modal, true);
                        end
                    end
                end
            end
        end
        
        function this=LabelBasics(contextProperty, ...
                contextDescription, nEvents, props, propsForLabels, ...
                dataSingular, dataPlural, subsetDescription)
            if nargin<8
                subsetDescription='leaf gate';
                if nargin<7
                    dataPlural='data';
                    if nargin<6
                        dataSingular='datum';
                        if nargin<5
                            propsForLabels=JavaProperties;
                            if nargin<4
                                props=JavaProperties;
                            end
                        end
                    end
                end
            end
            this.contextExcludeProperty=[LabelBasics.PROP_EXCLUDE contextProperty];
            this.contextDescription=contextDescription;
            this.contextHtml=['<font color="blue"><i>' ...
                contextDescription '</i></font>'];
            this.subsetDescription=subsetDescription;
            this.nEvents=nEvents;
            this.idPerEvent=zeros(this.nEvents, 1);
            
            this.propsForLabels=propsForLabels;
            this.props=props;
                        
            this.dataPlural=dataPlural;
            this.dataSingular=dataSingular;
            this.descriptions=JavaProperties;
        end
        
        function nEvents=addClass(this, id, ...
                sampleRows, name, count, color256)
            assert(isnumeric(id));
            assert(this.nEvents==length(sampleRows));
            if size(sampleRows,1)==1
                sampleRows=sampleRows';
            end
            if any(this.idPerEvent(:, end)~=0 & sampleRows)
                [~, cols]=size(this.idPerEvent);
                needNewOverFlow=true;
                if cols>1
                    for i=1:cols-1
                        wasSetAlready2=this.idPerEvent(:, i)~=0&sampleRows;
                        if  ~any(wasSetAlready2)
                            needNewOverFlow=false;
                            this.idPerEvent(sampleRows,i)=id;
                            break;
                        end
                    end
                end
                if needNewOverFlow
                    classIds=zeros(this.nEvents, 1);
                    classIds(sampleRows)=id;
                    this.idPerEvent=[this.idPerEvent classIds ];
                    this.idPerEvent(sampleRows,end)=id;
                end                
            else
                this.idPerEvent(sampleRows,end)=id;
            end
            nEvents=sum(sampleRows);
            this.doneGates=this.doneGates+1;
            id=num2str(id);
            this.descriptions.set(id, [name ' ' String.encodeK(count)]);
            this.propsForLabels.set(java.lang.String(id), name);
            key=this.ColorKey(id);
            if nargin>5
                this.propsForLabels.set(key, color256);
            else
                color256=this.props.get(key);
                if ~isempty(color256)
                    this.propsForLabels.set(key, color256);
                end
            end
        end
    end
end
