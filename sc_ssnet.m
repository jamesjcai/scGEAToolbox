function [A]=sc_ssnet(X,i,varargin)
% Construct single-sample network (ssnet) for given cell i
%
% see also: SC_GRN (single-sample GRN network)
% ref: https://github.com/WilfongGuo/Benchmark_control/tree/master/Benchmark_code
% https://doi.org/10.1371/journal.pcbi.1008962

p = inputParser;
defaultType = 'CSN';
validTypes = {'CSN','LIONESS','SSN','SPCC'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addRequired(p,'i',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,i,varargin{:})

   
switch upper(p.Results.type)
    case 'CSN'
        [A]=run.csnnet_ssnet(X,i);
        A=A{i};
    case 'LIONESS'             %  
        [A]=lioness_method(X,i);
    case 'SSN'      
        [A]=SSN(X(:,i),X);
    case 'SPCC'
        [A]=spcc_method(X,i);
end
if ~issparse(A)
    A=sparse(A);
end
end

function [CSN] = spcc_method(data,i)
%construct sample specific network

mean_a=mean(data,2);
std_a=std(data,0,2);
    v=data(:,i);
    x=(v-mean_a)./std_a;
    x(isnan(x))=0;
    x(abs(x)==inf)=0;
    
    cand1=(x*x');
    x=cand1;
    xx=triu(x);
    a_x=abs(xx(xx~=0));
    threshold=mean(a_x)+2*std(a_x);
    x(abs(x)<threshold)=0;
    x(x~=0)=1;
    CSN=x;




end


function [ index_R,p ] = SSN( sample,ref )
%function:construct the SSN
%   Input:
%         sample:calculated sample
%          ref:the reference samples
%   Output:
%         adjacency_matrix:the network structure
%a example
% sample=new_T(:,1);
% ref=new_N;

[R,P]=corrcoef(ref');
final_R0=R;
final_R0(isnan(final_R0))=0;

NEW_data=[ref sample];
[R1,P1]=corrcoef(NEW_data');
final_R1=R1;
final_R1(isnan(final_R1))=0;

index_R=final_R1-final_R0;
[m,n]=size(ref);
Z=index_R./((1-final_R0.^2)/(n-1));
Z(Z==inf)=max(max(Z));
Z(Z==-inf)=-max(max(Z));
Z(isnan(Z))=0;

clear NEW_data final_R1 final_R0 R0 R1 P P1
p=1-normcdf(abs(Z));

end


function [CSN] = lioness_method(data,i)
%construct sample specific network


PCC=corrcoef(data');
v0=data;

    v0(:,i)=[];
    PCC1=corrcoef(v0');
    x=size(data,2)*(PCC-PCC1)+PCC1;
    x(isnan(x))=0;
    
    xx=triu(x);
    a_x=abs(xx(xx~=0));
    threshold=mean(a_x)+2*std(a_x);
    
    x(abs(x)<threshold)=0;
    x(x~=0)=1;
  
    CSN=x;




end
