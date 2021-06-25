function [A]=sc_grn(X,varargin)
%SC_GRN construct single-cell gene regulatory network (scGRN)
%
% see also: SC_SSNET (single-sample network)

p = inputParser;
defaultType = 'pcnetpar';
validTypes = {'pcnet','pcnetpar','pcnetdenoised','genie3'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

   
switch p.Results.type
    case 'pcnet'
        [A]=sc_pcnet(X);
    case 'pcnetpar'             % parallel version 
        [A]=sc_pcnetpar(X);
    case 'pcnetdenoised'        % slow
        [A]=sc_pcnetdenoised(X);
    case 'genie3'               % slow
        [A]=run.GENIE3(X,[],true);
end
end
