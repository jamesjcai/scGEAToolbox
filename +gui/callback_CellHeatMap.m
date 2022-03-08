function callback_CellHeatMap(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);


% species=questdlg('Which species?','Select Species','Mouse','Human','Mouse');
% switch lower(species)
%     case 'human'
%         stag='hs';
%     case 'mouse'
%         stag='mm';
%     otherwise
%         return;
% end

f=figure;
[~,idx]=sort(sce.c);
X=sce.X(:,idx);

imagesc(log10(1+log10(1+X(1:1000,:))));
xlabel("Genes")
ylabel("Cells")

end
