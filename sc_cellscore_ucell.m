function [score] = sc_cellscore_ucell(X, genelist, tgsPos, tgsNeg)
%SC_CELLSCORE_UCELL   Cell-level gene signature scoring (UCell-style).
%  score = SC_CELLSCORE_UCELL(X, genelist, tgsPos, tgsNeg, opts)
%
%  X         : G × C expression matrix (genes × cells).
%  genelist  : G × 1 string/cell array of gene names.
%  tgsPos    : positive marker genes (string array)
%  tgsNeg    : negative marker genes (string array)
%
% Uses Mann–Whitney U normalization similar to UCell.
%
% see also: PKG.E_CELLSCORES, SC_CELLSCORE_ADMDL

% Compute cell scores from a list of feature genes
%
% tgsPos - positive features (positive target marker genes)
% tgsNeg - negative features (negative target marker genes)
%  opts.maxRank   : maximum rank threshold (default 1500)
%  opts.method    : 1 for rank-based, 2 for normalized expression


if nargin < 4 || isempty(tgsNeg)
    tgsNeg = ["IL2", "TNF"];
end
if nargin < 3
    tgsPos = ["CD44", "LY6C", "KLRG1", "CTLA", "ICOS", "LAG3"];
end

if ~any(matches(genelist, tgsPos, 'IgnoreCase', true))
    score = NaN(size(X, 2), 1);
    warning('No feature genes found in GENELIST.');
    return;
end

genelist = upper(genelist);
tgsPos = upper(tgsPos);
tgsNeg = upper(tgsNeg);

idx1 = matches(genelist, tgsPos, 'IgnoreCase', true);
idx2 = matches(genelist, tgsNeg, 'IgnoreCase', true);

n1 = sum(idx1);

if issparse(X)
    try
        X = full(X);
    catch
        warning('Could not convert sparse matrix to full. Proceeding with sparse.');
    end
end


methodid = 1;
switch methodid
    case 1
        % (UCELL-inspired Rank Method)
        % https://doi.org/10.1016/j.csbj.2021.06.043
        R = tiedrank(-X);
        R(R > 1500) = 1500 + 1;
        u = sum(R(idx1, :)) - (n1 * (n1 - 1)) / 2;
        score = 1 - u / (n1 * 1500);
        score = score(:);
    case 2
        % (Rank-Sum Based Scoring)
        X = sc_norm(X);
        disp('Library-size normalization...done.')
        X = log1p(X);
        disp('Log(x+1) transformation...done.')
        data_avg = mean(X, 2);
        idx1 = matches(genelist, tgsPos, 'IgnoreCase', true);
        x = X(idx1, :);
        [~, ~, stats0] = ranksum(mean(x, 2), data_avg);
        score = zeros(size(X, 2), 1);
        for k = 1:size(X, 2)
            [~, ~, stats] = ranksum(x(:, k), data_avg);
            score(k) = stats.ranksum - stats0.ranksum;
        end
        %score(score<0)=0;
        %score=log(exp(score)+1);
        score = (score - min(score)) ./ (max(score) - min(score));
end

end

%{


function STATS=mwwtest(x1,x2)
% Mann-Whitney-Wilcoxon non parametric test for two unpaired groups.
% This file execute the non parametric Mann-Whitney-Wilcoxon test to evaluate the
% difference between unpaired samples. If the number of combinations is less than
% 20000, the algorithm calculates the exact ranks distribution; else it
% uses a normal distribution approximation. The result is not different from
% RANKSUM MatLab function, but there are more output informations.
% There is an alternative formulation of this test that yields a statistic
% commonly denoted by U. Also the U statistic is computed.
%
% Syntax: 	STATS=MWWTEST(X1,X2)
%
%     Inputs:
%           X1 and X2 - data vectors.
%     Outputs:
%           - T and U values and p-value when exact ranks distribution is used.
%           - T and U values, mean, standard deviation, Z value, and p-value when
%           normal distribution is used.
%        If STATS nargout was specified the results will be stored in the STATS
%        struct.
%
%      Example:
%
%         X1=[181 183 170 173 174 179 172 175 178 176 158 179 180 172 177];
%
%         X2=[168 165 163 175 176 166 163 174 175 173 179 180 176 167 176];
%
%           Calling on Matlab the function: mwwtest(X1,X2)
%
%           Answer is:
%
% MANN-WHITNEY-WILCOXON TEST
%
%                        Group_1    Group_2
%                        _______    _______
%
%     Numerosity          15         15
%     Sum_of_Rank_W      270        195
%     Mean_Rank           18         13
%     Test_variable_U     75        150
%
% Sample size is large enough to use the normal distribution approximation
%
%     Mean       SD        Z       p_value_one_tail    p_value_two_tails
%     _____    ______    ______    ________________    _________________
%
%     112.5    24.047    1.5386    0.061947            0.12389
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2009). MWWTEST: Mann-Whitney-Wilcoxon non parametric test for two unpaired samples.
% http://www.mathworks.com/matlabcentral/fileexchange/25830

%Input Error handling
p = inputParser;
addRequired(p,'x1',@(x) validateattributes(x,{'numeric'},{'row','real','finite','nonnan','nonempty'}));
addRequired(p,'x2',@(x) validateattributes(x,{'numeric'},{'row','real','finite','nonnan','nonempty'}));
parse(p,x1,x2);

%set the basic parameter
n1=length(x1); n2=length(x2); NP=n1*n2; N=n1+n2; N1=N+1; k=min([n1 n2]);

[A,B]=tiedrank([x1(:); x2(:)]); %compute the ranks and the ties
R1=A(1:n1); R2=A(n1+1:end);
T1=sum(R1); T2=sum(R2);
U1=NP+(n1*(n1+1))/2-T1; U2=NP-U1;
disp('MANN-WHITNEY-WILCOXON TEST')
disp(' ')
disp(table([n1;T1;T1/n1;U1],[n2;T2;T2/n2;U2],...
    'VariableNames',{'Group_1' 'Group_2'},...
    'RowNames',{'Numerosity' 'Sum_of_Rank_W' 'Mean_Rank' 'Test_variable_U'}))
if nargout
    STATS.n=[n1 n2];
    STATS.W=[T1 T2];
    STATS.mr=[T1/n1 T2/n2];
    STATS.U=[U1 U2];
end
if round(exp(gammaln(N1)-gammaln(k+1)-gammaln(N1-k))) > 20000
    mU=NP/2;
    if B==0
        sU=realsqrt(NP*N1/12);
    else
        sU=realsqrt((NP/(N^2-N))*((N^3-N-2*B)/12));
    end
    Z1=(abs(U1-mU)-0.5)/sU;
    p=1-normcdf(Z1); %p-value
    disp('Sample size is large enough to use the normal distribution approximation')
    disp(' ')
    disp(table(mU,sU,Z1,p,2*p,'VariableNames',{'Mean' 'SD' 'Z' 'p_value_one_tail' 'p_value_two_tails'}))
    if nargout
        STATS.method='Normal approximation';
        STATS.mU=mU;
        STATS.sU=sU;
        STATS.Z=Z1;
        STATS.p=[p 2*p];
    end
else
    disp('Sample size is small enough to use the exact Mann-Whitney-Wilcoxon distribution')
    disp(' ')
    if n1<=n2
        w=T1;
    else
        w=T2;
    end
    pdf=sum(nchoosek(A,k),2);
    P = [sum(pdf<=w) sum(pdf>=w)]./length(pdf);
    p = min(P);
    disp(table(w,p,2*p,'VariableNames',{'W' 'p_value_one_tail' 'p_value_two_tails'}))
    if nargout
        STATS.method='Exact distribution';
        STATS.T=w;
        STATS.p=[p 2*p];
    end
end
disp(' ')
end


%}