function A = sc_grn(X, type)
    %SC_GRN Construct single‑cell gene regulatory network (scGRN)
    %
    %   A = sc_grn(X) uses default type 'pcnetpar'.
    %   A = sc_grn(X, type) where type is one of:
    %       'pcnet', 'pcnetpar', 'pcnetdenoised', 'genie3'
    
    arguments
        X {mustBeNumeric}
        type (1,1) string = "pcnetpar"
        % or char, depending on preference
        % Uncomment next line if you want to allow char inputs too
        % type {mustBeMember(type,["pcnet" "pcnetpar" "pcnetdenoised" "genie3"])}
    end
    
    type = lower(string(type));
    validTypes = ["pcnet" "pcnetpar" "pcnetdenoised" "genie3"];
    if ~ismember(type, validTypes)
        error("sc_grn:InvalidType", ...
              "Type must be one of: %s", strjoin(validTypes, ", "));
    end
    
    switch type
        case "pcnet"
            A = sc_pcnet(X);
        case "pcnetpar"
            A = sc_pcnetpar(X);
        case "pcnetdenoised"
            A = ten.sc_pcnetdenoised(X);
        case "genie3"
            % GENIE3 can be slow—consider future support for parfeval if needed
            A = run.GENIE3(X, [], true);
    end
end


%{
function [A] = sc_grn(X, varargin)
%SC_GRN construct single-cell gene regulatory network (scGRN)
%
% see also: SC_SSNET (single-sample network)

p = inputParser;
defaultType = 'pcnetpar';
validTypes = {'pcnet', 'pcnetpar', 'pcnetdenoised', 'genie3'};
checkType = @(x) any(validatestring(x, validTypes));

addRequired(p, 'X', @isnumeric);
addOptional(p, 'type', defaultType, checkType)
parse(p, X, varargin{:})
    
    switch p.Results.type
        case 'pcnet'
            [A] = sc_pcnet(X);
        case 'pcnetpar' % parallel version
            [A] = sc_pcnetpar(X);
        case 'pcnetdenoised' % slow
            [A] = ten.sc_pcnetdenoised(X);
        case 'genie3' % slow
            [A] = run.GENIE3(X, [], true);
    end
end
%}