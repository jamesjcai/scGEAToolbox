function [K,usehvgs] = i_gethvgnum(sce)
K=[];
usehvgs=false;
                answer2 = questdlg(sprintf('Use highly variable genes (HVGs, n=2000) or use all genes (n=%d)?', sce.NumGenes), ...
                    '', '2000 HVGs 🐇', 'All Genes 🐢', 'Other...', '2000 HVGs 🐇');
                switch answer2
                    case 'All Genes 🐢'
                        usehvgs = false;
                        K = sce.NumGenes;
                    case '2000 HVGs 🐇'
                        usehvgs = true;
                        K = 2000;
                    case 'Other...'
                        K = gui.i_inputnumk(min([3000, sce.NumGenes]), ...
                            100, sce.NumGenes);
                        if isempty(K), return; end
                        usehvgs = true;
                    otherwise
                        return;
                end

