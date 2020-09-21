function [ class,type ] = DRSCAN(x,threshold,Eps,Eta,Matrix)
% run ReCon-DBSCAN based on a dissmilarity matrix
%   Input  : X  : data set (m,n); m-objects, n-variables 
%            threshold :  density-ratio threshold for DBSCAN
%            Eps :  Eps neigbourhood for ReCon-DBSCAN
%            Eta :  Eta neigbourhood for ReCon-DBSCAN    
%            Matrix :  dissmilarity matrix
%   Output : class: cluster labels (m by 1)
%            type :  label types (1: core; 2: boundary; 3: noise)


[m,n]=size(x);

M=Matrix-Eps;
b=(M<=0);
DenEps=sum(b,2);

M=Matrix-Eta;
b=(M<=0);
DenEta=sum(b,2);

%Ratio=DenEps./DenEta *(Eta/Eps)^n; % find core points
Ratio=DenEps./DenEta; % find core points

type=zeros(1,m);
type(Ratio>=threshold)=1;

%% linking core points
x=[[1:m]' x];
touched=zeros(m,1);
class=zeros(1,m);
no=1;
for i=1:m
    if touched(i)==0;
        ob=x(i,:);
        D=Matrix(i,:);
        ind=find(D<=Eps);
        if type(i)==1 && length(ind)>1
            class(ind)=ones(length(ind),1)*max(no);
            while ~isempty(ind)
                obi=ind(1);
                touched(ind(1))=1;
                class(obi)=no;
                ind(1)=[];
                
                if type(obi)==1 
                    D=Matrix(obi,:);
                    i1=find(D<=Eps);
                    class(i1)=no;
                    
                    for i=1:length(i1)
                        if touched(i1(i))==0
                            touched(i1(i))=1;
                            ind=[ind i1(i)];
                            class(i1(i))=no;
                        end
                    end
                    
                end
            end
            no=no+1;
        end
    end
end

i1=find(class==0);
class(i1)=-1;
type(i1)=-1;
end


