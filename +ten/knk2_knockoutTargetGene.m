function [T]=knk2_knockoutTargetGene(A0,targetgene,genelist,dosort,lambdav)

%see also: ten.i_knk
    import ten.*
    if nargin<5, lambdav=0; end
    if nargin<4, dosort=true; end
    idx=find(genelist==targetgene,1);
    if isempty(idx)
        error("TARGETGENE is not found in GENELIST");
    end
    [T]=i_knk(A0,idx,genelist,dosort,lambdav);
end
