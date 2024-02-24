function [glist] = i_selectngenes(sce, predefinedlist, parentfig)

% This function uses i_selectngenes
if nargin < 2, predefinedlist = []; end
if nargin < 3, parentfig = []; end

% internal function used by callback_BuildGeneNetwork
glist = [];
if isa(sce, 'SingleCellExperiment')
    gsorted = natsort(sce.g);
elseif isstring(sce)
    genelist = sce;
    gsorted = natsort(genelist);
end


if ~isempty(predefinedlist)
    predefinedlist = gsorted(matches(gsorted, predefinedlist, 'IgnoreCase', true));
end

answer = questdlg('Select genes from list or paste gene names?', ...
    'Select/Paste Genes', 'Select', 'Paste', 'Cancel', 'Select');
switch answer
    case 'Cancel'
        return;
    case 'Paste'
        n = length(gsorted);
        if isempty(predefinedlist)
            tg = gui.i_inputgenelist(gsorted(randperm(n, min([20, length(gsorted)]))));
        else
            tg = gui.i_inputgenelist(predefinedlist);
        end

        if length(tg) >= 1
            [y, ix] = ismember(upper(tg), upper(gsorted));
            % i=i(y);
            %glist=tg(y);
            glist = gsorted(ix(y));
            a = length(glist) - length(tg);
            if a ~= 0
                if a == 1
                    warndlg(sprintf('%d gene is not found.', a));
                elseif a > 1
                    warndlg(sprintf('%d genes are not found.', a));
                end
            end
            %             if length(glist)<2
            %                 warndlg('Need at leaset 2 genes');
            %                 return;
            %             end
        else
            %             warndlg('Need at least 2 genes');
            return;
        end
    case 'Select'

        if isa(sce, 'SingleCellExperiment')
            [gsorted] = gui.i_sortgenenames(sce);
            if isempty(gsorted), return; end
        end

        [idx] = gui.i_selmultidlg(gsorted, predefinedlist, parentfig);
        if isempty(idx), return; end
        if idx == 0, return; end
        %if length(idx)<2
        %    warndlg('Need at least 2 genes');
        %    return;
        %else
        glist = gsorted(idx);
        %g='Dhfr, Lmbr1, Reck, Rnf168, Rpl26, Snrnp27, Tmem160'
        %g=["Tcf7","Lef1","Bcl6","Ctla4","Lag3","Pdcd1"];
        %end
    otherwise
        return;
end

end