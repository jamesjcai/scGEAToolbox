function callback_ShowGeneExpr(src, ~)

if ismcc || isdeployed
    makePPTCompilable();
    % https://www.mathworks.com/help/rptgen/ug/compile-a-presentation-program.html
end
import mlreportgen.ppt.*;

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
[axx, bxx] = view(findall(FigureHandle,'type','axes'));

% answer = 'Multiple';
% switch answer
%     case 'Single'
%         [gsorted] = gui.i_sortgenenames(sce);
%         if isempty(gsorted), return; end
%         [indx, tf] = listdlg('PromptString', {'Select a gene:', '', ''}, ...
%             'SelectionMode', 'single', 'ListString', gsorted);
%         if tf == 1
%             [methodid] = gui.i_pickscatterstem('Scatter+Stem');
%             % methodid=2;      case 'Scatter'
%             % 'Stem'           methodid=1;
%             if isempty(methodid), return; end
%             for k = 1:length(indx)
%                 gui.i_cascadefig(sce, gsorted(indx(k)), axx, bxx, k, methodid);
%             end
%         end
%     case 'Multiple'
        [glist] = gui.i_selectngenes(sce);
        if isempty(glist), return; end
        sc_uitabgrpfig(sce,glist);
end



%{

            [y, i] = ismember(upper(glist), upper(sce.g));
            if ~all(y), error('Unspecific running error.'); end
            glist = sce.g(i);

            %[~,i]=ismember(gsorted(idx),sce.g);
            x = sum(sce.X(i, :), 1);
            if length(i) == 1
                g = sce.g(i);
            elseif length(i) > 1
                answer2 = questdlg('Intersection (AND) or Union (OR)', ...
                    '', 'Individually', ...
                    'Intersection (AND)', ...
                    'Union (OR)', ...
                    'Individually');
                switch answer2
                    case 'Union (OR)'
                        g = sprintf("%s | ", glist);
                    case 'Intersection (AND)'
                        g = sprintf("%s & ", glist);
                        ix = sum(sce.X(i, :) > 0, 1) == length(i);
                        if ~any(ix)
                            % helpdlg('No cells express all the selected genes.', '');
                            return;
                        end
                        x = x .* ix;
                    case 'Individually'


                        [methodid] = gui.i_pickscatterstem('Scatter+Stem');
                        if isempty(methodid), return; end

                        answer = questdlg('Output to PowerPoint?');
                        switch answer
                            case 'Yes'
                                needpptx = true;
                            case 'No'
                                needpptx = false;
                            otherwise
                                return;
                        end

                        images = cell(length(glist), 1);
                        warning off
                        for k = 1:length(glist)
                            f = gui.i_cascadefig(sce, glist(k), axx, bxx, k, methodid);
                            % i_showcascade(sce,gsorted(idx(k)),axx,bxx,k);
                            if needpptx
                                img1 = [tempname, '.png'];
                                images{k} = img1;                                
                                saveas(f, img1);
                            end
                        end
                        warning on
                        if needpptx, gui.i_save2pptx(images); end
                        return;
                    otherwise
                        return;
                end % end of AND / OR / Individual
                g = extractBefore(g, strlength(g)-2);
            end
            f = figure('visible', 'off');
            [h1] = sc_scattermarker(x, g, sce.s, g, 5);
            title(g);
            view(h1, axx, bxx);

            movegui(f, 'center');
            set(f, 'visible', 'on');

%}
