function [k,usehvgs] = i_gethvgnum(sce, parentfig)

if nargin<2, parentfig = []; end

k=[];
usehvgs=false;
    answer = gui.myQuestdlg(parentfig, ...
        sprintf(['Use highly variable genes (HVGs, n=2000) ' ...
        'or use all genes (n=%d)?'], sce.NumGenes), ...
        '', {'2000 HVGs 🐇', 'All Genes 🐢', 'Other...'}, '2000 HVGs 🐇');
    switch answer
        case 'All Genes 🐢'
            usehvgs = false;
            k = sce.NumGenes;
        case '2000 HVGs 🐇'
            usehvgs = true;
            k = 2000;
        case 'Other...'
            k = gui.i_inputnumk(min([3000, sce.NumGenes]), ...
                100, sce.NumGenes);
            if isempty(k), return; end
            usehvgs = true;
        otherwise
            return;
    end
end

