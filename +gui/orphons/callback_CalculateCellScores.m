function callback_CalculateCellScores(src, ~, sce)
% UNUSED function will be removed.

if nargin < 3
    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
end

actiontype = questdlg('Select a predefined score or define a new score?', ...
    '', 'Select Predefined Score', ...
    'Define New Score', 'Select Predefined Score');
switch actiontype
    case 'Select Predefined Score'
        [~, T] = pkg.e_cellscores([], [], 0);

        listitems = T.ScoreType;
        [indx, tf] = listdlg('PromptString', ...
            {'Select Class'}, ...
            'SelectionMode', 'single', ...
            'ListString', listitems, ...
            'ListSize', [220, 300]);
        if ~tf == 1, return; end


        fw = gui.gui_waitbar;
        [cs, ~, posg] = pkg.e_cellscores(sce.X, sce.g, indx);
        ttxt = T.ScoreType(indx);
        gui.gui_waitbar(fw);

    case 'Define New Score'
        ttxt = 'Customized Score';
        % Pd1=pdcd1 tim3=HAVCR2, tcf1=HNF1A  https://www.nature.com/articles/s41577-019-0221-9
        % posgcandidates=["PDCD1","HNF1A","HAVCR2","KLRG1","CD44","LY6C","CTLA","ICOS","LAG3"];
        %posgcandidates=sce.g(randi(length(sce.g),10,1));
        [posg] = gui.i_selectngenes(sce.g);
        if isempty(posg)
            helpdlg('No feature genes selected.', '')
            return;
        end
        %         fw=gui.gui_waitbar;
        %         %a=sprintf('%s+',posg);
        %         %a=a(1:min([length(a),50]));
        %         %ttxt=sprintf('%s\n%s',ttxt,a);
        %         posg
        %         cs=sc_cellscore_admdl(sce.X,sce.g,posg);
        %         gui.gui_waitbar(fw);

        [cs] = gui.e_cellscore(sce, posg);
        ttxt = {ttxt};

    otherwise
        return;
end


if ~isempty(cs)
    gui.i_stemscatterfig(sce, cs, posg, ...
        matlab.lang.makeValidName(ttxt{1}));
end


answer2 = questdlg(sprintf('Score has been computed.\nCompare the score across cell classes?'), 'Continue?');
switch answer2
    case 'Yes'

    otherwise
        return;
end
[thisc, clabel] = gui.i_select1class(sce);
if isempty(thisc) % || numel(unique(thisc))==1
    errordlg('Undefined');
    return;
end

gui.i_violinplot(y, thisc, ttxt);
xlabel(clabel);

%figure('WindowStyle','modal');
%pkg.i_violinplot_groupordered(cs,thisc);
%ylabel(strrep(ttxt,'_','\_'))


end
