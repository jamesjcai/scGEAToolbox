function [T, A0] = sctenifoldknk(X, genelist, kogene, varargin)
% T=sctenifoldknk(X,genelist,"Foxp3");
%
% X is a gene x cell matrix from the wild-type
%
% CHECK THE RESULT AGAINST A NULL BEFORE USING IT. T's p-values come from
% TEN.I_DR, which fits a chi-square to drdist/mean(drdist). That asks whether
% a gene moved more than the average gene moved in this run - not whether it
% moved because of KOGENE. On a network with strong hub structure the same
% genes top the ranking whichever row is zeroed, and the enrichment computed
% from them describes the network rather than the knockout.
%
% TEN.KNKNULLCONTROL settles it by knocking out random genes from the same A0
% and rescoring against that background. It costs minutes, needs the A0
% returned here, and reports a verdict:
%
%     [T, A0] = ten.sctenifoldknk(X, genelist, "Foxp3");
%     [Tnull, S] = ten.knknullcontrol(A0, "Foxp3", genelist);
%     disp(S.Verdict)
%
% see also: TEN.KNKNULLCONTROL, TEN.I_KNK, TEN.I_DR, TEN.SCTENIFOLDNET
import ten.*

if nargin < 3
    error(sprintf('USAGE: T=sctenifoldknk(X,genelist,kogene);\n       T=sctenifoldnet_m(X0,X1,genelist,''qqplot'',true);'));
end
if isscalar(kogene) && isnumeric(kogene)
    idx = kogene;
else
    idx = find(genelist == kogene, 1);
    if isempty(idx)
        error("KOGENE should be a member of GENELIST.");
    end
end
p = inputParser;
addOptional(p, 'sorttable', false, @islogical);
addOptional(p, 'smplmethod', "bootstrap", @(x) (isstring(x) | ischar(x)) & ismember(lower(string(x)), ["jackknife", "bootstrap"]));
addOptional(p, 'tdmethod', "CP", @(x) (isstring(x) | ischar(x)) & ismember(upper(string(x)), ["CP", "TUCKER"]));
addOptional(p, 'nsubsmpl', 10, @(x) fix(x) == x & x > 0);
addOptional(p, 'csubsmpl', 500, @(x) fix(x) == x & x > 0);
addOptional(p, 'savegrn', false, @islogical);
parse(p, varargin{:});
dosort = p.Results.sorttable;
tdmethod = p.Results.tdmethod;
nsubsmpl = p.Results.nsubsmpl;
csubsmpl = p.Results.csubsmpl;
smplmethod = p.Results.smplmethod;
savegrn = p.Results.savegrn;

switch upper(tdmethod)
    case "CP"
        tdmethod = 1;
    case "TUCKER"
        tdmethod = 2;
end
switch lower(smplmethod)
    case "jackknife"
        usebootstrp = false;
    case "bootstrap"
        usebootstrp = true;
end

if size(X, 1) ~= length(genelist)
    error('Length of genelist should be the same as the number of rows of X0 or X1.');
end

% Add the Tensor Toolbox from the saved preference if it is not already on
% the path, then error only if it genuinely is not installed.
%
% This used to error outright, which made the function fail in any MATLAB
% that had not happened to call ten.sctenifoldnet first - that one does the
% addpath itself, so an interactive session picked the path up as a side
% effect and a fresh `matlab -batch` did not. A long unattended run would
% then die at its first contrast. ten.check_tensor_toolbox is the toolbox's
% own helper for exactly this and is what ten.sctenifoldnet's own preamble
% amounts to.
ten.check_tensor_toolbox;
    %    if exist('ten.i_td1.m','file')~=2
    %        error('Need i_td1.m in the scTendifoldNet https://github.com/cailab-tamu/scTenifoldNet/tree/master/MATLAB');
    %    end
    %    if exist('sc_pcnet.m','file')~=2
    %        error('Need sc_pcnet.m in the scGEAToolbox https://github.com/jamesjcai/scGEAToolbox');
    %    end

X = sc_norm(X, "type", "libsize");
X = log1p(X);

[XM] = i_nc(X, nsubsmpl, 3, csubsmpl, usebootstrp);
[A0] = i_td1(XM, tdmethod);
    if savegrn
        tstr = matlab.lang.makeValidName(strrep(sprintf("GRN_created_on_% s", datetime)," ", "_at_"));
        save(tstr, 'A0', 'genelist', '-v7.3');
        fprintf('\nConstructed gene regulatory network (GRN) is saved in A0_%s.mat\n', tstr);
    end
    %     A1=A0;
    %     A1(idx,:)=0;
    %     [aln0,aln1]=i_ma(A0,A1);
    %     T=i_dr(aln0,aln1,genelist,dosort);
T = ten.i_knk(A0, idx, genelist, dosort);
end
