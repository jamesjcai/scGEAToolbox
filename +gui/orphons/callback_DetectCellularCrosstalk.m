function callback_DetectCellularCrosstalk(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[thisc, ~] = gui.i_select1class(sce, false);
if isempty(thisc), return; end

[c, cL] = grp2idx(thisc);

[idx] = gui.i_selmultidlg(cL);
if isempty(idx), return; end
if numel(idx) < 2
    warndlg('Need at least 2 cell types');
    return;
end
selected = ismember(c, idx);
fw = gui.gui_waitbar;
sce = sce.selectcells(selected);


[OUT, T] = run.talklr(sce.X, sce.g, cL(c(selected)));
gui.gui_waitbar(fw);


n = length(OUT.ligandok);
if n == 0
    warndlg('Not detected.');
    return;
end

gui.i_exporttable(T, true, 'T');

%     if ~(ismcc || isdeployed)

%     labels = {'Save OUT to variable named:'};
%     vars = {'OUT'};
%     values = {OUT};

%         [f,ft]=export2wsdlg(labels, vars, values);
%         waitfor(f);
%         if ft
%             disp('Run >> gui.i_crosstalkgraph(OUT,k,sce); to plot crosstalk graph for ligand-receptor pair k.')
%         end
%     end

answer = questdlg('Plot top ligand-receptor expression?');

switch answer
    case 'Yes'
        %[numfig]=gui.i_inputnumg(n,'the number of ligand-receptor pairs');
        [numfig] = gui.i_inputnumk(10, 1, n, 'the number of ligand-receptor pairs');

        if isempty(numfig) || isnan(numfig), return; end
        if isnumeric(numfig) && numfig > 0 && numfig <= n
            F = cell(numfig, 1);
            %listitems=cell(numfig,1);
            for k = 1:numfig
                F{k} = gui.i_crosstalkgraph(OUT, k, sce);
                %listitems{k}=sprintf('%s -> %s (KL = %.2f)',...
                %    OUT.ligandok(k), OUT.receptorok(k),...
                %    OUT.KL(k));
                drawnow;
            end
            gui.i_export2pptx(F);
            %gui.i_export2pptx(F,string(listitems));
        end
end


% answer=questdlg('Interactive exploration?');
% switch answer
%     case 'Yes'
%         listitems=cell(n,1);
%         for k=1:n
%             listitems{k}=sprintf('%s -> %s (KL = %.2f)',...
%                 OUT.ligandok(k), OUT.receptorok(k),...
%                 OUT.KL(k));
%         end
%         i_displyres(listitems);
% end
%
%     function i_displyres(listitems)
%         [indx2,tf2] = listdlg('PromptString',...
%             {'Select ligand-receptor pairs to plot'},...
%              'SelectionMode','single','ListString',listitems,...
%              'ListSize',[210,300]);
%          if tf2==1
%                 kk=indx2;
%                 gui.i_crosstalkgraph(OUT,kk,sce);
%                 i_displyres(listitems);
%          elseif tf2==0
%
%          end
%     end


end