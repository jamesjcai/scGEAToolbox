function [glist] = i_selectngenes(sce,predefinedlist)

if nargin<2, predefinedlist=[]; end

% internal function used by callback_BuildGeneNetwork
glist=[];
if isa(sce,'SingleCellExperiment')
    gsorted=sort(sce.g);
elseif isstring(sce)
    genelist=sce;
    gsorted=sort(genelist);
end

if ~isempty(predefinedlist)
    predefinedlist=gsorted(matches(gsorted,predefinedlist,'IgnoreCase',true));
end

answer = questdlg('Select genes from list or paste gene names?',...
	'Select/Paste Genes','Select','Paste','Cancel','Select');
switch answer
    case 'Cancel'
        return;
    case 'Paste'
        n=length(gsorted);
        if isempty(predefinedlist)
            tg=gui.i_inputgenelist(gsorted(randperm(n,20)));
        else
            tg=gui.i_inputgenelist(predefinedlist);
        end
        if length(tg)>=2
            
            [y]=ismember(upper(tg),upper(gsorted));
            % i=i(y);
            glist=tg(y);
%             if length(glist)<2
%                 warndlg('Need at leaset 2 genes');
%                 return;
%             end
        else
%             warndlg('Need at least 2 genes');
            return;
        end
    case 'Select'
        [idx]=gui.i_selmultidlg(gsorted,predefinedlist);
        if isempty(idx), return; end
        %if length(idx)<2
        %    warndlg('Need at least 2 genes');
        %    return;
        %else
            glist=gsorted(idx);
            %g='Dhfr, Lmbr1, Reck, Rnf168, Rpl26, Snrnp27, Tmem160'
            %g=["Tcf7","Lef1","Bcl6","Ctla4","Lag3","Pdcd1"];
        %end
end

end