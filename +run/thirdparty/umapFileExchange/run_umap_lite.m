function [reduction]=run_umap_lite(varargin)

[args, argued]=UmapUtil.Initialize(varargin{:});
args=UmapUtil.CheckArgs(args, argued);
%clusterIdentifiers=[];
%globals=BasicMap.Global;
%homeFolder=globals.appFolder;
% try
%     props=fullfile(homeFolder, 'globals.mat');
%     globals.load(props);
% catch
% end

% reductions=Map;
% reductions.load(fullfile(homeFolder, 'reductions.mat'));
% args.reduction=reductions.newId;

reduction=[];
umap=[];


%csv_file_or_data=args.csv_file_or_data;
save_template_file = args.save_template_file;
beQuiet=true;
extras=UMAP_extra_results;

% args=UmapUtil.RelocateExamples(args);

csv_file_or_data=args.csv_file_or_data;

inData=csv_file_or_data;
parameter_names=args.parameter_names;

template_file=args.template_file;
[nRows, nCols]=size(inData);
 
% sCols=num2str(nCols);
% if nRows*nCols<15
%     if length(inData)==1
%         if inData>1&&inData<=19
%             if askYesOrNo(['Run example #' num2str(inData) '??'])
%                 run_examples(inData)
%                 return
%             end
%         end
%     end
%     if isempty(inData)
%     else
%         msgError(sprintf(...
%             'Too little data: %d row(s) X  %d col(s)??!', ...
%             nRows, nCols));
%     end
%     reduction=[];
%     umap=[];
%     clusterIdentifiers=[];
%     extras=[];
%     return
% end
% if strcmpi(args.label_column, 'end')
%     args.label_column=nCols;
% end
% if args.label_column>0 
%     if  args.label_column>nCols
%         % msg(Html.WrapHr(['The input data has ' sCols ' columns ...<br>'...
%         %     'THUS the label_column must be >=1 and <= ' sCols]));
%         % assert(args.label_column>0 && args.label_column<=nCols, [...
%         %     'label_column must be >=1 and <= ' sCols]);
%     end
%     labelCols=1;
%     if args.matchingTestLabels
%         testSetLabels=inData(:, args.label_column);
%     end
%     % if ~isempty(template_file) 
%     %     if args.label_column<=length(parameter_names)
%     %         parameter_names(args.label_column)=[];
%     %     end
%     % else
%     %     labels=inData(:, args.label_column);        
%     % end
%     inData(:, args.label_column)=[];
% else
%     labelCols=0;
% end


umap = UMAP;

% isSupervising=isprop(umap, 'supervisors') && ~isempty(umap.supervisors);
% 
% if isSupervising
%     if args.n_components==2
%         progressMatchType=0;
%     else
%         limit=args.match_3D_limit;
%         if limit<nRows
%             progressMatchType=-1;
%         else
%             progressMatchType=3;
%         end
%     end
% else
%     if any(args.match_supervisors>0)
%         if argued.contains('match_supervisors')
%             if ~isequal(1, args.match_supervisors)
%                 warning('match_supervisors only affects supervised template reduction');
%                 args.match_supervisors=1;
%             end
%         end
%     end
% end

UmapUtil.SetArgsTemplateCanOverride(umap, args, argued, parameter_names);
umap.n_epochs=args.n_epochs;
umap.nn_descent_min_rows=args.nn_descent_min_rows;
umap.nn_descent_min_cols=args.nn_descent_min_cols;
umap.nn_descent_max_neighbors=args.nn_descent_max_neighbors;
umap.nn_descent_transform_queue_size=args.nn_descent_transform_queue_size;

        
%method=umap.setMethod(args.method);
umap.verbose=~beQuiet;
umap.randomize=args.randomize;


%labelMap=[];
nParams=length(parameter_names);

% good=nParams==0||nParams==nCols || (args.label_column>0 &&...
%     (nParams==nCols-1 || nParams==nCols));


% if ~good
%     if args.label_column>0
%         preAmble=sprintf(['# data columns=%d, # parameter_names=%d '...
%             'since label_column=%d <br># parameter_names must be '...
%             '%d or %d '],  nCols, nParams, args.label_column, ...
%             nCols, nCols-1);
%     else
%         preAmble=sprintf(['# of data columns(%s) must equal '...
%             '# of parameter_names(%d)'], sCols, nParams);
%     end
%     % msg(Html.WrapHr(preAmble));
%     % assert(nParams==0||nParams==nCols || (args.label_column>0 &&...
%     %     (nParams==nCols-1 || nParams==nCols)), preAmble);    
% end

% if ~isempty(newSubsetIdxs)
%     hasLabels=true;
%     labelCols=0;
%     [labels, labelMap]=resupervise(umap, inData, newSubsetIdxs);
%     nLabels=length(unique(labels));
% elseif args.label_column>0 && isempty(template_file) && ~args.matchingUmap
%     hasLabels=true;
%     labelCols=1;
%     good=args.label_column>0 && args.label_column<=nCols;
%     if ~good
%         msg(Html.WrapHr(['The input data has ' sCols ' columns ...<br>'...
%             'THUS the label_column must be >=1 and <= ' sCols]));
%         assert(args.label_column>0 && args.label_column<=nCols, [...
%             'label_column must be >=1 and <= ' sCols]);
% 
%     end
%     nLabels=length(unique(labels));
%     if nLabels > .5*nRows
%         preAmble='%d is a LOT of distinct labels for a %dx%d matrix!';
%         msg(['WARNING:  ' sprintf(preAmble, nLabels, nRows, nCols)]);
%         warning(preAmble, nLabels, nRows, nCols);
%     end
%     if args.label_column<=nParams
%         parameter_names(args.label_column)=[];
%     end
%     umap.dimNames=parameter_names;
%     nLabels=length(unique(labels));
% else

%hasLabels=false;

%    nLabels=0;
% end


% if args.label_column>0 
%     if isSupervising && ~isempty(args.label_file) && isempty(testSetLabels)
%         if args.label_column==0
%             warning(['label_file has no effect when reducing '...
%                 'with supervised template without test set labels']);
%         elseif isempty(find(args.match_scenarios==1, 1)) .... 
%                 && isempty(find(args.match_scenarios>2,1))
%             warning(['test set labels only needed for umap'...
%                 ' supervised templates IF specifying '...
%                 'match_scenarios 1 3 or 4 ']);
%         end
%     else
%         if ~isempty(testSetLabels)
%             [labelMap, halt]=getLabelMap(testSetLabels);
%         else
%             [labelMap, halt]=getLabelMap(labels);
%             if hasLabels
%                 ColorsByName.Override(labelMap, args.color_file, beQuiet);
%             end
%         end
%         if halt
%             delete(fig);
%             globals.save;
%             return;
%         end
%         if ~isempty(labelMap) && ~args.buildLabelMap ...
%                 && isSupervising && ~isempty(testSetLabels)
%             testSetLabels=LabelBasics.RelabelIfNeedBe(testSetLabels, ...
%                 sprv.labelMap, labelMap);
%         end
% 
%     end
% end


% nanRows=any(isnan(inData'));
% badRows=sum(nanRows);
% if badRows>0
%         warning(['Data matrix has ' ...
%                 String.Pluralize2('row', badRows) ...
%                 'with NAN values!']);
%         inData=inData(~nanRows,:);
%     if any(isnan(inData(:)))
%         showMsg(Html.WrapHr(['Sorry...<br>can not proceed<br>'...
%             '<br>NAN values exist... SIGH!']));
%         globals.save;
%         return;
%     end
% end

% args.hiD=nCols-labelCols;
% info=[String.encodeInteger(nRows) 'x' String.encodeInteger(nCols-labelCols)];
% if ischar(csv_file_or_data)
%     [~, fileName]=fileparts(csv_file_or_data);
%     info=['UMAP on ' String.ToLaTex(fileName) ', ' info];
% else
%     info=['[UMAP on ' info];
% end
% if args.python
%     info=[info ', Python'];
% else
%     info=[info ', ' method];
% end
% pause(.01);
% info2=['(optimize\_layout method=' method ')'];




%READY TO START REDUCING HiD to LoD!!
% if isSupervising
%     args.reductionType=UMAP.REDUCTION_SUPERVISED_TEMPLATE;
% else    
    %if isempty(template_file)
        % if hasLabels
        %     % args.reductionType=UMAP.REDUCTION_SUPERVISED;
        % else
        %     args.reductionType=UMAP.REDUCTION_BASIC;
        % end
    %else
    %    args.reductionType=UMAP.REDUCTION_TEMPLATE;
    %end
%end
%extras.args=args;
%strReduction=[UmapUtil.GetReductionLongText(args.reductionType) ' reduction'];
%reportProgress(['Running ' strReduction ', v' UMAP.VERSION], true);
%Map.SetStruct(reductions, args.reduction, args);

% reductions.save;

%    if ~args.python
            reduction = umap.fit_transform(inData);
            % assignin("base","xumap",umap);
%    end
    % if ~isempty(umap.supervisors)
    %     sprv=umap.supervisors;
    %     sprv.description=args.description;
    % end

% if ~isempty(reduction)
%     testBasicReduction;
%     reportProgress(['Finished ' strReduction]);
% else
%     msg('Parameter reduction was cancelled or not done');
%     if exist('pu', 'var') && isa(pu, 'PopUp')
%         pu.stop;
%         pu.dlg.dispose;
%     end
% end

% if (nargout>1 || ~isempty(save_template_file)) && ~isempty(reduction)
% 
%     umap.prepareForTemplate;
% 
% end    

% if nargout>2 || ~strcmpi(args.cluster_output, 'none')
%     if nargout<3 && ~strcmpi(args.cluster_output, 'graphic')
%         warning('No clusterIdentifiers output argument');
%     elseif ~strcmpi(args.cluster_output, 'ignore')
%         clusterIdentifiers=doClusters(reduction);
%         if isempty(clusterIdentifiers)
%             dispNoDbScan;
%         end
%     end
% end
% if ~beGraphic
%     if nargout<4 %don't delete figures if expecting them in extras
%         extras.closeMatchFigs;
%         extras.closeTreeFigs;
%     end
% end
%globals.save;
 
end
