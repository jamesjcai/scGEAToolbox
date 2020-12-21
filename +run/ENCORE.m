function ENCORE(X,genelist)
if nargin<1, X=rand(4,5); end
if nargin<2, genelist=string((1:4)'); end

[isok,msg]=commoncheck_R('R_ENCORE');
if ~isok, error(msg); end

writematrix(X,'input.txt');
writematrix(genelist,'genelist.txt');

% cd(oldpth);
end