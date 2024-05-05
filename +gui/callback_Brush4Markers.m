function callback_Brush4Markers(src, event)
if exist('lasso', 'file')
    % disp('using LASSO')
    i_Brush4MarkersLASSO(src, event);
    return;
end
end


function i_Brush4MarkersLASSO(src, ~, sce)
FigureHandle = src.Parent.Parent;

if nargin < 3
    sce = guidata(FigureHandle);
end

%    FigureHandle=src.Parent.Parent;
%    sce=guidata(FigureHandle);
%   assert(isequal(FigureHandle.Children,...
%         FigureHandle.findobj('type','Axes')))

% axesh=FigureHandle.Children(1)
axesh = FigureHandle.findobj('type', 'Axes');
[axx, bxx] = view(axesh);

% [axx, bxx] = view(findall(FigureHandle,'type','axes'));

assert(isequal(axesh.findobj('type', 'Scatter'), ...
    FigureHandle.findobj('type', 'Scatter')))
%axesh.Children(1)
%isequal(axesh.findobj('type','Scatter'),axesh.Children(2))
h = axesh.findobj('type', 'Scatter');
ptsSelected = logical(h.BrushData.');


if ~any(ptsSelected)
    warndlg("No cells are brushed/selected.",'','modal');
    answer=questdlg('Select cells by a grouping variable?','');
    if ~strcmp(answer,'Yes'), return; end
    [ptsSelected] = gui.i_select1classcells(sce, false);    
    if isempty(ptsSelected), return; end
    if all(ptsSelected)
        warndlg("All cells are in the same group.",'');
        return;
    end
else
    % assignin('base', 'ptsSelected', ptsSelected);
    [ptsSelected, letdoit] = gui.i_expandbrushed(ptsSelected, sce);
    if ~letdoit, return; end
end


[numfig] = gui.i_inputnumg(500);
if isempty(numfig), return; end
fw = gui.gui_waitbar;
y = double(ptsSelected);
sce.c = 1 + ptsSelected;
X = sce.X';
try
    if issparse(X), X = full(X); end
    [B] = lasso(X, y, 'DFmax', numfig*3, 'MaxIter', 1e3);
catch ME
    gui.gui_waitbar(fw);
    errordlg(ME.message);
    rethrow(ME);
end

[~, ix] = min(abs(sum(B > 0)-numfig));
b = B(:, ix);
idx = b > 0;
gui.gui_waitbar(fw);

if ~any(idx)
    warndlg('No marker found','')
    return;
end


    markerlist = sce.g(idx);
    [~, jx] = sort(b(idx), 'descend');
    markerlist = markerlist(jx);

    fprintf('%d marker genes: ', length(markerlist));
    fprintf('%s ', markerlist)
    fprintf('\n')

%    gui.i_exporttable(table(markerlist), true, ...
%        'Tmarkerlist','MarkerListTable');

    % 'Tviolindata','ViolinPlotTable'
    % 'Tmarkerlist','MarkerListTable'

%    [answer] = questdlg('Plot expression of markers?');
%    if isempty(answer), return; end
%    switch answer
%        case 'Yes'
            % [methodid] = gui.i_pickscatterstem('Scatter');
            % if isempty(methodid), return; end
            % F = cell(length(markerlist), 1);
            % for kk = 1:length(markerlist)
            %     F{kk} = gui.i_cascadefig(sce, markerlist(end-(kk - 1)), ...
            %         ax, bx, kk, methodid);
            % end
            % gui.i_export2pptx(F, flipud(markerlist(:)));
%        fw=gui.gui_waitbar;
        gui.sc_uitabgrpfig_expplot(sce, markerlist, FigureHandle, [axx, bxx]);
        gui.gui_waitbar(fw);    
%    end
%end

%     pause(2);
%     export2wsdlg({'Save marker list to variable named:'},...
%             {'g_markerlist'},{markerlist});
end


%{


FigureHandle=src.Parent.Parent;
sce=guidata(FigureHandle);
%     assert(isequal(FigureHandle.Children, FigureHandle.findobj('type','Axes')))
%
%     axesh=FigureHandle.Children(1);
axesh=FigureHandle.findobj('type','Axes');
[ax,bx]=view(axesh);
assert(isequal(axesh.findobj('type','Scatter'),...
    FigureHandle.findobj('type','Scatter')))
h=axesh.Children(1);
ptsSelected = logical(h.BrushData.');

if ~any(ptsSelected)
    warndlg("No cells are selected.");
    return;
end
assignin('base','ptsSelected',ptsSelected);

[c,cL]=grp2idx(sce.c);
if isscalar(unique(c))
    methodtag=1;
else
    answer = questdlg('Select brushed cell group?');
    if strcmp(answer,'Yes')
        if isscalar(unique(c(ptsSelected)))
            methodtag=2;
        else
            errordlg('More than one group of brushed cells');
            return;
        end
    elseif strcmp(answer,'No')
        methodtag=1;
    else
        return;
    end
end
[numfig]=gui.i_inputnumg;
if isempty(numfig), return; end
fw=gui.gui_waitbar;

switch methodtag
    case 1
        [markerlist]=sc_pickmarkers(sce.X,sce.g,1+ptsSelected,2);
        sce.c=1+ptsSelected;
        markerlist=markerlist{2};
        % disp('xxx')
    case 2
        ptsSelected=c==unique(c(ptsSelected));
        % h.BrushData=double(ptsSelected);
        %[markerlist]=sc_pickmarkers(sce.X,sce.g,c,unique(c(ptsSelected)));
        [markerlist]=sc_pickmarkers(sce.X,sce.g,1+ptsSelected,2);
        sce.c=1+ptsSelected;
        markerlist=markerlist{2};
end
gui.gui_waitbar(fw);
% assignin('base','A',A);
[numfig]=gui.i_inputnumg;
fw=gui.gui_waitbar;
htmlfilename=cL{unique(c(ptsSelected))};
pkg.i_markergeneshtml(sce,markerlist,numfig,...
    [ax bx],htmlfilename,ptsSelected);
gui.gui_waitbar(fw);
%     pause(2);
%     export2wsdlg({'Save marker list to variable named:'},...
%             {'g_markerlist'},{markerlist});
end
%}
