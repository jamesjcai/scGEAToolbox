function [ranks,tail] = calc_ranks(data,group,method,ifabs,iftrank)
% Function for calculating gene ranking metric for selected method
% When external function is used, it must take 2 input arguments 
% (data matrix and sample labels) and return ranking values.
% Author:
% Michal.Marczyk@polsl.pl

n = length(group);
if size(data,2)~= n
    data = data';
end
m = size(data,1);
gr0 = group==0; n0 = sum(gr0);
gr1 = group==1; n1 = sum(gr1);

if isa(method, 'function_handle')
    m_fun = method;
    method = 'external';
end

switch method
    case 'S2N'
        ranks = (mean(data(:,gr1),2) - mean(data(:,gr0),2))./(std(data(:,gr1),[],2) + std(data(:,gr0),[],2));
        tail = 'both';
    case 'ttest'
        denom = sqrt(var(data(:,gr1),[],2)/n1 + var(data(:,gr0),[],2)/n0);
        ranks =  (mean(data(:,gr1),2) - mean(data(:,gr0),2))./(std(data(:,gr1),[],2) + std(data(:,gr0),[],2))./denom;
        tail = 'both';
    case 'ratio'
        ranks = mean(2.^data(:,gr1),2)./mean(2.^data(:,gr0),2);
        tail = 'both';
    case 'diff'
        ranks = mean(data(:,gr1),2) - mean(data(:,gr0),2);
        tail = 'both';
    case 'log2_ratio'
        ranks = log2(mean(2.^data(:,gr1),2)./mean(2.^data(:,gr0),2));
        tail = 'both';
    case 'SoR'
        ranks = zeros(1,m);
        for a=1:m
            rank_tmp = tiedrank(data(a,:));
            ranks(a) = sum(rank_tmp(gr1));
        end
        tail = 'right';
    case 'BWS'
        ranks = zeros(1,m);
        ind0 = 1:n0;
        ind1 = 1:n1;
        for a=1:m
            rank_tmp = tiedrank(data(a,:));
            R1 = sort(rank_tmp(gr0));
            R2 = sort(rank_tmp(gr1));
            B0 = sum((R1 - ((n1+n0)/n0 * ind0)).^2./((ind0./(n0+1)) .* (1-(ind0/(n0+1))) * ((n1*(n1+n0))/n0)));
            B1 = sum((R2 - ((n1+n0)/n1 * ind1)).^2./((ind1./(n1+1)) .* (1-(ind1/(n1+1))) * ((n0*(n1+n0))/n1)));
            ranks(a) = (B0/n0 + B1/n1)/2;
        end
        tail = 'right';
    case 'ReliefF'
        [~,ranks] = relieff(data',group,10);
        tail = 'right';
    case 'WAD'
        ranks = WAD(data,group);
        tail = 'both';
    case 'FCROS'
        ranks = FCROS(data,group);
        tail = 'both';
    case 'MWT'
        ranks = mwt(data,group);
        tail = 'both';
    case 'MSD'
        ranks = MSD(data,group);
        tail = 'right';
    case 'external'
        ranks = m_fun(data,group);
        tail = 'both';
    otherwise
        error('Wrong gene ranking method selected.')
end
ranks = ranks(:);

if ifabs
    ranks = abs(ranks);
    tail = 'right';
end

if iftrank
    ranks = tiedrank(ranks);
    tail = 'right';
end

function WAD_stat = WAD(data,group)
group_un = unique(group);   %find unique gruoups

%calculate average difference
dm_1 = mean(data(:,group==group_un(1)),2);
dm_2 = mean(data(:,group==group_un(2)),2);
AD = dm_2 - dm_1;

%calculate weights and WAD
dm_all = (dm_1 + dm_2)/2;
w = (dm_all - min(dm_all))/range(dm_all);
WAD_stat = w .* AD;

function res = FCROS(data,group)
d1 = data(:,group==1);
n_d1 = sum(group==1);

d0 = data(:,group==0);
[n,n_d0] = size(d0);

rank=zeros(n,n_d1*n_d0);
it=1;

for i=1:n_d1
    for j=1:n_d0
        c = d1(:,i)./d0(:,j);
        rank(:,it) = tiedrank(c)/n;
        it = it + 1;
    end
end
res = trimmean(rank,10,2);

function W = mwt(x,grp)
grp_un = unique(grp); 
n1 = sum(eq(grp_un(1),grp)); n2 = sum(eq(grp_un(2),grp));   %group sizes

if ~eq(size(x,2),n1+n2)
    x = x';
end

if ~eq(size(x,2),n1+n2)
    error('Wrong number of samples.')
end

%means and standard deviations of gene expressions
mi1 = mean(x(:,eq(grp_un(1),grp)),2); mi2 = mean(x(:,eq(grp_un(2),grp)),2);
var1 = var(x(:,eq(grp_un(1),grp)),[],2); var2 = var(x(:,eq(grp_un(2),grp)),[],2);

d = n1 + n2 - 2;    %df of the pooled standard error

%pooled and unpooled standard error
sep2 = (((n1-1)*var1 + (n2-1)*var2)/d)*(1/n1 + 1/n2);
seu2 = var1/n1 + var2/n2;

%df of the unpooled standard error
df = (seu2.*seu2)./(((1/(n1-1))*(var1/n1).*(var1/n1)) + ((1/(n2-1))*(var2/n2).*(var2/n2)));

%posterior prob. that the group variances are equal
F = var1./var2;
pF = 2*min(fcdf(F,n1-1,n2-1),fcdf(1./F,n2-1,n1-1));
warning('off','bioinfo:mafdr:PoorEstimatedPI0Value')
w = mafdr(pF);
warning('on','bioinfo:mafdr:PoorEstimatedPI0Value')

%weighted standard error
sew2 = w.*sep2 + (1-w).*seu2;
dw = w.*d + (1-w).*df;

%hyperparameters estimated from the data
[s02,d0] = est_hyper(log(sew2),nanmean(dw));

%moderated standard error and modified Welch statistic
sem = (d0.*s02 + dw.*sew2)./(d0 + dw);    %moderated standard error
W = (mi1 - mi2)./sqrt(sem);    %moderated Welch statistic

function [s02,d0] = est_hyper(z,D)

z(isnan(z)) = 0;
fun1 = @(d0,z,D) var(z) - psi(1,D/2) - psi(1,d0/2);
limit = fun1(100,z,D);
if limit<0
    d0 = 100;
else
    fun2 = @(d0) fun1(d0,z,D);
    d0 = fzero(fun2,[0.01,100]);
end

s02 = exp(mean(z) - psi(0,D/2) + psi(0,d0/2) - log(d0/D));

function stat = MSD(x,grp)
grp_un = unique(grp); 
n1 = sum(eq(grp_un(1),grp)); n2 = sum(eq(grp_un(2),grp));   %group sizes

%means and standard deviations of gene expressions
mi1 = mean(x(:,eq(grp_un(1),grp)),2); mi2 = mean(x(:,eq(grp_un(2),grp)),2);
var1 = var(x(:,eq(grp_un(1),grp)),[],2); var2 = var(x(:,eq(grp_un(2),grp)),[],2);

%calculate moderated T
d = n1 + n2 - 2;    %df of the pooled standard error
sp2 = (((n1-1)*var1 + (n2-1)*var2)/d);   %pooled standard error
[s02,d0] = est_hyper(log(sp2),nanmean(d));     %hyperparameters estimated from the data
se = sqrt((1/n1 + 1/n2)*(d0.*s02 + d.*sp2)./(d0 + d));    %moderated standard error

%calculate log fold change and 95% CI
lFC = mi2-mi1;
ind = lFC > 0;

stat = zeros(size(lFC));
stat(ind) = lFC(ind) - tinv(.975,d+d0)*se(ind);
stat(~ind) = -lFC(~ind) - tinv(.975,d+d0)*se(~ind);

% tmp = 68;
% tstat = lFC./se;
% [lFC(tmp),tstat(tmp),stat(tmp),se(tmp),d0,lFC(tmp)-tinv(.975,d+d0)*se(tmp),lFC(tmp)+tinv(.975,d+d0)*se(tmp)]