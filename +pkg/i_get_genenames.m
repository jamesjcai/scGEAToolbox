function [T]=i_get_genenames(g)
%Get gene names

    pw=fileparts(mfilename('fullpath'));
    dbfile=fullfile(pw,'..','resources','HGNCBiomart.txt');
    Tright=readtable(dbfile);
    g=upper(g(:));
    if ~iscell(g)
        ApprovedSymbol=cellstr(g);
    else
        ApprovedSymbol=g;
    end
    Tleft=table(ApprovedSymbol);
    T = join(Tleft,Tright);
    
