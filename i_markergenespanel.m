function i_markergenespanel(X,genelist,s,markerlist,numfig)

if nargin<5, numfig=5; end
for kkk=1:numfig
    figure;
    for kk=1:min([9,length(markerlist)])
        subplot(3,3,kk)
        sc_scattermarker(X,genelist,s,...
            markerlist(kk+9*(kkk-1)),3,5,false);
    end
end
