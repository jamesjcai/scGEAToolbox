function callback_BuildGeneNetwork(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    gsorted=sort(sce.g);
answer = questdlg('Select genes to construct regulatory network?',...
	'Build scGRN');
switch answer
    case 'Yes'
        [idx]=gui.gui_selmultidlg(gsorted);
        if isempty(idx)
            helpdlg('No gene selected.');
        else            
            [~,i]=ismember(gsorted(idx),sce.g);
            sprintf("%s ",gsorted(idx))
            %g='Dhfr, Lmbr1, Reck, Rnf168, Rpl26, Snrnp27, Tmem160'
            %g=["Tcf7","Lef1","Bcl6","Ctla4","Lag3","Pdcd1"];
            %[~,i]=ismember(g,sce.g);
            fw=gui.gui_waitbar;    
            x=sce.X(i,:);
            A=sc_pcnet(x);
%             B=e_transf(A,0.6);            
%             G=digraph(B,gsorted(idx),'omitselfloops');
            gui.gui_waitbar(fw);            
            sc_grnview(A,gsorted(idx));
        end
end
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