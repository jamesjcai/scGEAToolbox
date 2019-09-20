[X,genelist,celltypeid]=sc_readtsvfile('example_data/yan.csv');
rng(235); showlegend=true;
% rng(111); showlegend=false;
% rng(113); showlegend=true;

s=sc_tsne(X,2);
c1=sc_sc3(X,6);
c2=run_simlr(X,6);
c3=run_soptsc(X,'k',6);


% i_myscatter(s,grp2idx(celltypelist))
fh=figure; 
subplot(2,2,1)
gscatter(s(:,1),s(:,2),celltypelist)
%legend('Location','northeastoutside')
if showlegend, legend('Location','northwest'); else, legend off; end
title('Cell type (''Ground Turth'')')

subplot(2,2,2)
gscatter(s(:,1),s(:,2),c1)
%legend('Location','southwestoutside')
if showlegend, legend('Location','northwest'); else, legend off; end
title('SC3')

subplot(2,2,3)
gscatter(s(:,1),s(:,2),c2)
%legend('Location','southwestoutside')
if showlegend, legend('Location','northwest'); else, legend off; end
title('SIMLR')

subplot(2,2,4)
gscatter(s(:,1),s(:,2),c3)
%legend('Location','southwestoutside')
if showlegend, legend('Location','northwest'); else, legend off; end
title('SoptSC')

fh.Position=[fh.Position(1) fh.Position(2)-100 fh.Position(3)+100 fh.Position(4)+100];
