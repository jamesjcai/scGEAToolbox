function callback_BuildGeneNetwork(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    gsorted=sort(sce.g);
answer = questdlg('Paste or select genes?',...
	'Build scGRN','Paste','Select','Cancel','Paste');
switch answer
    case 'Cancel'
        return;
    case 'Paste'
        n=length(sce.g);
        tg=gui.gui_inputgenelist(sce.g(randi(n,10,1)));        
        if length(tg)>=2
            [y,i]=ismember(tg,sce.g);
            i=i(y);
            glist=tg(y);
            if length(glist)<2
                warndlg('Need at leaset 2 genes');
                return;
            end
        else
            warndlg('Need at least 2 genes');
            return;
        end
    case 'Select'
        [idx]=gui.gui_selmultidlg(gsorted);
        if isempty(idx), return; end
        if length(idx)<2
            warndlg('Need at least 2 genes');
            return;
        else
            glist=gsorted(idx);
            %g='Dhfr, Lmbr1, Reck, Rnf168, Rpl26, Snrnp27, Tmem160'
            %g=["Tcf7","Lef1","Bcl6","Ctla4","Lag3","Pdcd1"];
        end
end
            [y,i]=ismember(glist,sce.g);
            if ~all(y), error('xxx'); end
            fprintf("%s\n",glist)
            fw=gui.gui_waitbar;    
            x=sce.X(i,:);
            A=sc_pcnet(x);
%             B=e_transf(A,0.6);            
%             G=digraph(B,gsorted(idx),'omitselfloops');
            gui.gui_waitbar(fw);            
            sc_grnview(A,glist);
end


function A=e_transf(A,q)
% A - adjacency matrix
if nargin<2, q=0.95; end
dim=size(A);
if numel(dim)==2
    a=max(abs(A(:)));
    if a>0
        A=A./a;
        A=A.*(abs(A)>quantile(abs(A(:)),q));        
    end
elseif numel(dim)==3
    for k=1:dim(3)
        A(:,:,k)=e_transf(A(:,:,k),q);
    end
end
end