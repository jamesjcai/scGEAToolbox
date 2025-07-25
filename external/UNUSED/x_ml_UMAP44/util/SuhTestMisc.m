classdef SuhTestMisc
%   This class is a "grab bag" of miscellaneous static test functions

%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
    
%#ok<*NOPRT,*ASGLU>
    methods(Static)

        function FlowJoLabelGates(labels, data, columnNames, gates)
            FlowJoTree.CreateLabelGates('PHATE',...
                    'demoIcon.gif',data, labels, [], columnNames, gates);
        end

        function clues=Cluster(data, ax)
            [~, clues, dbm]=Density.GetClusters(data, 'medium', min(data), max(data));
            dbm.drawBorders(ax);
        end

        function name=Name(btn)
            app=BasicMap.Global;
            name=char(btn.getText);
            idx=String.IndexOf(name, app.supStart);
            if idx>0
                name=name(1:idx-1);
            end
            name=Html.Remove(name);
            name=strrep(name, '&bull;', '');
        end

        function  [data, names, labelPropsFile, gates, leaves, sampleOffsets]...
                =packageSubsets(this, justData, askUser, ids, purpose)
            if nargin<5
                purpose='UMAP';
                if nargin<4
                    ids=this.getSelectedIds;
                    if nargin<3
                        askUser=true;
                    end
                end
            end
            nIds=length(ids);
            fullNameMap=[];
            if nIds>1 && askUser
                [yes, cancelled]=askYesOrNo(...
                    ['<html>Run ' purpose ' on <br>' ...
                        String.Pluralize2('gate selection', nIds) ...
                        '?<hr></html>']);
                if cancelled
                    return;
                end
                if ~yes
                    ids(2:end)=[];
                    nIds=1;
                else
                    fullNameMap=Map;
                end
            end
            gates={};
            if nIds==1
                [data, names, labelPropsFile, ~, gates{1}, leaves]...
                    =this.packageSubset(ids{1}, false, ...
                    this.pipelineFcsParameters, false, justData, fullNameMap);
            else
                fcsNames=this.getCommonColumnNames(ids, askUser, this.pipelineFcsParameters);
                [data, names, labelPropsFile, ~, gates{1}, leaves, props]...
                    =this.packageSubset(ids{1}, false, ...
                    fcsNames, false, justData, fullNameMap);
                resave=false;
                for i=2:nIds
                    [data2, ~, ~, ~, gate2, leaves2, props2]...
                        =this.packageSubset(ids{i}, askUser, ...
                        fcsNames, false, justData, fullNameMap);
                    if ~isempty(data2)
                        gates{end+1}=gate2;
                        data=[data;data2];
                        leaves=[leaves leaves2];
                        if ~justData
                            keys=props.keys;
                            nKeys=length(keys);
                            for j=1:nKeys
                                if ~props.containsKey(keys{j})
                                    resave=true;
                                    props.set(key, props2.get(keys{j}))
                                end
                            end
                        end
                    end
                end
                if resave
                    props.save;
                end
            end
            sampleOffsets=Map;
            sampleOffset=1;
            gateOffset=1;
            for i=1:nIds
                gate=gates{i};
                nEvents=length(gate.getSampleRows);
                offsets=[sampleOffset, gateOffset, nEvents, gate.count];
                sampleOffsets.set(gate.id, offsets);
                sampleOffset=sampleOffset+nEvents;
                gateOffset=gateOffset+gate.count;
            end
        end
        
        function ok=fncCheck(jd, finalAnsw)
            ok=true;
            fprintf('jd="%s" and finalAnsw="%s"\n', jd.getTitle, finalAnsw);
        end

        function concludeUmap(this, ~, ~)
            fjw=this.gml;
        end

        function [X, Y]=newUmapXY(this)
            if nargin<3
                ask=true;
            end
            originalName=this.getName;
            gml=this.gml;
            next=1;
            if ~ask
                gateName=originalName;
            end
            [~, names]=gml.getDerived(this.gater.sampleNum);
            nSibs=length(names);
            X=[];
            Y=[];
            while true
                if ask
                    gateName=inputDlg(struct('where', 'north+', ...
                        'msg', 'Enter a UMAP name/identifer  ...'), ...
                        'New derived parameter ...', originalName);
                    if isempty(gateName)
                        return;
                    end
                end
                found=false;
                for i=1:nSibs
                    X=['UMAP_1_' gateName]; Y=['UMAP_2_' gateName];
                    if strcmpi(X, names{i}) || strcmpi(Y, names{i})
                        if ask
                            if ~askYesOrNo(struct('icon', 'warning.png',...
                                    'msg', Html.WrapHr(sprintf([...
                                    '"<b>%s</b>"<br>is already used ' ...
                                    'in this sample for UMAP<br><br>Retry?'], ...
                                    [X '/' Y]))), 'DUPLICATE!', 'north+')
                                return;
                            end
                        else
                            next=next+1;
                            gateName=[originalName ' #' num2str(next)];
                        end
                        found=true;
                        break;
                    end
                end
                if ~found
                    break;
                end
            end
            
        end

        function [toNodes, fromNodes, names, columnIndexes, files, types, minRanges, maxRanges]...
                =copyDervived(gml, fromSampleNum, toSampleNum, name)
            if nargin<4
                name='UMAP*';
                if nargin<3
                    toSampleNum=1;
                    if nargin<2
                        fromSampleNum=2;
                    end
                end
            end
            [fromNodes, names, columnIndexes, files, types, minRanges, ...
                maxRanges, gains]=gml.getDerived(fromSampleNum, name);
            toNodes={};
            for i=1:length(fromNodes)
                toNodes{end+1}=gml.addDerivedLinearCsv(toSampleNum, ...
                    columnIndexes(i), files{i}, names{i}, maxRanges(i), ...
                    minRanges(i), gains(i));
            end
        end
        function [nodes, names, columnIndexes, files, types, minRanges, maxRanges]...
                =derivedTest(gml, sampleNum)

            [nodes, names, columnIndexes, files, types, minRanges, maxRanges]=gml.getDerived(sampleNum, 'UMAP_2_8GQ0')
            [nodes, names, ~,~,~,minRanges,maxRanges]=gml.getDerived(sampleNum, 'UMAP*')
            [nodes, names, columnIndexes, files, types, minRanges, maxRanges]=gml.getDerived(sampleNum)
        end

        function refreshChildren(this)
            if ~isempty(this.gater) && ~isempty(this.gater.tree)
                this.gater.gml.resyncChildren(this.id);
                this.gater.tree.suhTree.ensureChildUiNodesExist(this.id);
            end
        end
        function [popNode, gateNode, id]=createSubGate(gate, roiType, ...
                roiPosition, name, dims, count, scalers)
            if nargin<7
                scalers=[];
            end
            gml=gate.gml;
            doc=gml.doc;
            [yes, subPopNode]=gml.hasSubPopulations(gate.population);
            if ~yes
                subPopNode=doc.createElement(gml.xmlSubpop);
                gate.population.appendChild(subPopNode);
            end
            
            popNode=doc.createElement(gml.xmlPop);
            popNode.setAttribute('name', name);
            popNode.setAttribute('count', num2str(count));
            popNode.setAttribute('annotation', '');
            popNode.setAttribute('owningGroup', '');
            popNode.setAttribute('expanded', '0');
            popNode.setAttribute('sortPriority', '10');
            subPopNode.appendChild(popNode);
            gml.addParentGraph(popNode, gate.population, dims);
            gateNode=doc.createElement(gml.xmlGate);
            pid=gml.ParseId(gate.id);
            id=gml.newGateId;
            setIds(gateNode);
            popNode.appendChild(gateNode);
            [gatingNodeName, gatingType]=gml.getGatingNode(roiType);
            gatingNode=doc.createElement(gatingNodeName);
            setIds(gatingNode);
            setDim(dims{1});
            setDim(dims{2});
            gateNode.appendChild(gatingNode);
            if ~isempty(scalers)
                scalerX=scalers{1};
                scalerY=scalers{2};
            else
                scalerX=gate.fcs.scalers.get(dims{1});
                scalerY=gate.fcs.scalers.get(dims{2});
            end
            if strcmp(gatingType, gml.xmlPolygon)
                [R, C]=size(roiPosition);
                assert(C==2);
                for i=1:R
                    setCoordinates(roiPosition(i, 1), roiPosition(i,2));
                end
            end
            gml.save(true);
            
            function setIds(node)
                node.setAttribute(gml.xmlPid, pid);
                node.setAttribute(gml.xmlId, id);
            end

            function setDim(dim)
                dim1=doc.createElement(gml.xmlDim);
                gatingNode.appendChild(dim1);
                dim2=doc.createElement(gml.nameFcsDim);
                dim2.setAttribute(gml.xmlTypeName, dim)
                dim1.appendChild(dim2);
            end

            function setCoordinates(X, Y)
                vertex=doc.createElement(gml.nameVertex);
                gatingNode.appendChild(vertex);
                coord=doc.createElement(gml.nameCoord);
                coord.setAttribute(gml.nameValue, ...
                    num2str(scalerX.inverse(X)))
                vertex.appendChild(coord);
                coord=doc.createElement(gml.nameCoord);
                coord.setAttribute(gml.nameValue, ...
                    num2str(scalerY.inverse(Y)))
                vertex.appendChild(coord);
            end
        end

        function changes=clearMissingGates(this, idSet)
            props=this.props;
            ids=props.keys;
            N=length(ids);
            changes=0;
            for i=1:N
                id=ids{i};
                if startsWith(id, FlowJoWsp.TYPE_GATE)
                    strs=strtrim(props.get(id));
                    toks=java.lang.String(id).split('[.:]');
                    if isempty(strs) 
                        continue;
                    end
                    if ~idSet.contains(toks(2))
                        props.remove(id);
                        changes=changes+1;
                        continue;
                    end
                elseif startsWith(id, FlowJoWsp.TYPE_SAMPLE)
                    strs=strtrim(props.get(id));
                else
                    continue;
                end
                ids2=strsplit(strs, ' ');
                N2=length(ids2);
                for j=1:N2
                    id2=ids2{j};
                    if startsWith(id2, FlowJoWsp.TYPE_GATE)
                        toks=java.lang.String(id2).split('[.:]');
                        if ~idSet.contains(toks(2))
                            props.remove(id{j});
                            changes=changes+1;
                            break;
                        end
                    else
                        disp('huh');
                    end
                end
            end
            
        end
        
        function changes=clearGates(this)
            props=this.props;
            ids=props.keys;
            N=length(ids);
            changes=0;
            for i=1:N
                id=ids{i};
                if startsWith(id, FlowJoWsp.TYPE_GATE)
                        changes=changes+1;
                        props.remove(id);
                elseif startsWith(id, FlowJoWsp.TYPE_SAMPLE)
                    if endsWith(id, '.children')
                    changes=changes+1;
                    props.remove(id);
                    end
                end
            end
            
        end
    end
end