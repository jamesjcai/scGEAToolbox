function [thisc, clabel] = i_select1class(sce, allowsingle, ...
    promptstr, ...
    prefersel, parentfig)

if nargin < 5, parentfig = []; end
if nargin < 4, prefersel = ''; end
if nargin < 3 || isempty(promptstr)
    promptstr = 'Select grouping variable:';
end
if nargin < 2 || isempty(allowsingle), allowsingle = true; end
if ~isempty(parentfig)
    figure(parentfig);
    cleanupObj = onCleanup(@() figure(parentfig));
end
thisc = [];
clabel = '';

[listitems, a, b] = i_classlistitems(sce, allowsingle);

% listitems={'Current Class (C)','Cluster ID','Batch ID',...
%            'Cell Type','Cell Cycle Phase'};

y = false;
if ~isempty(prefersel)
    [y, idx]=ismember(prefersel,listitems);
end

if gui.i_isuifig(parentfig)
    % allowmulti = false, explicitly. myListdlg defaults it to TRUE, which
    % used to give this dialog a multi-select listbox that the code below
    % cannot honour: listitems{indx2} with a vector indx2 quietly yields
    % only the first value, so extra picks were accepted and then dropped
    % without a word. The listdlg fallback right below has always said
    % 'single'; this is the uifigure branch agreeing with it, and with the
    % name of this function. To grouping by several variables at once, use
    % gui.i_selectnclass.
    [indx2, tf2] = gui.myListdlg(parentfig, listitems, ...
        promptstr, prefersel, false);
else
    if y
        [indx2, tf2] = listdlg('PromptString', ...
            {promptstr}, ...
            'SelectionMode', 'single', 'ListString', listitems, ...
            'ListSize', [220 300], 'InitialValue', idx, 'Name', ' ');
    else
        [indx2, tf2] = listdlg('PromptString', ...
            {promptstr}, ...
            'SelectionMode', 'single', 'ListString', listitems, ...
            'ListSize', [220 300], 'Name', ' ');
    end
end


if tf2 == 1
    clabel = listitems{indx2};
    switch clabel
        case 'Current Class (C)'
            thisc = sce.c;
        case 'Cluster ID' % cluster id
            thisc = sce.c_cluster_id;
        case 'Batch ID' % batch id
            thisc = sce.c_batch_id;
        case 'Cell Type' % cell type
            thisc = sce.c_cell_type_tx;
        case 'Cell Cycle Phase' % cell cycle
            thisc = sce.c_cell_cycle_tx;
        case 'Workspace Variable...'
            thisc = i_pickvariable;
    end
else
    figure(parentfig);
end


function [c] = i_pickvariable
        c = [];
        %     a=evalin('base','whos');
        %     b=struct2cell(a);
        %     v=false(length(a),1);
        %     for k=1:length(a)
        %         if max(a(k).size)==sce.NumCells && min(a(k).size)==1
        %             v(k)=true;
        %         end
        %     end
        %     if any(v)
        % valididx=ismember(b(4,:),'double');
        % a=a(valididx);
        if gui.i_isuifig(parentfig)
            % Same reason as above, and a nastier failure if left multi:
            % a(indx).name with a vector indx expands to a comma-separated
            % list, turning line 'evalin(''base'', name)' into the
            % three-argument evalin(context, expr, catch_expr) form - which
            % silently evaluates the second pick only when the first throws.
            [indx, tf] = gui.myListdlg(parentfig, b(1, :), ...
                'Select variable:', [], false);
        else
            [indx, tf] = listdlg('PromptString', {'Select variable:'}, ...
                'liststring', b(1, :), 'SelectionMode', 'single');
        end
        if tf == 1
            c = evalin('base', a(indx).name);
        end
        %    end
end
    % if isempty(thisc)
    %     gui.myErrordlg(parentfig, 'Undefined');
    %     return;
    % end
    % if numel(unique(thisc))==1
    %     gui.myWarndlg(parentfig, "Cannot compare with an unique group");
    %     return;
    % end
end
