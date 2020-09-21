function [class,type]=Mdbscan(x,k,Eps,Matrix)
% run DBSCAN based on a dissmilarity matrix
%   Input  : X  : data set (m,n); m-objects, n-variables 
%            k :  MinPts for DBSCAN
%            Eps :  eps neigbourhood for DBSCAN
%            Matrix :  dissmilarity matrix
%   Output : class: cluster labels (m by 1)
%            type :  label types (1: core; 2: boundary; 3: noise)

[m,n]=size(x);

if nargin<3 || isempty(Eps)
    [Eps]=epsilon(x,k);
end

x=[[1:m]' x];
[m,n]=size(x);
type=zeros(1,m);
no=1;
touched=zeros(m,1);
class=zeros(1,m);
for i=1:m
    if touched(i)==0;
        ob=x(i,:);
        %%
        %D=dist(ob(2:n),x(:,2:n));  %1*m matrix
        D=Matrix(i,:);
        %%
        ind=find(D<=Eps);
        
        if length(ind)>1 & length(ind)<k
            type(i)=0;
            class(i)=0;
        end
        if length(ind)==1
            type(i)=-1;
            class(i)=-1;
            touched(i)=1;
        end
        
        if length(ind)>=k;
            type(i)=1;
            class(ind)=ones(length(ind),1)*max(no);
            
            while ~isempty(ind)
                ob=x(ind(1),:);
                touched(ind(1))=1;
                
                %%
                %D=dist(ob(2:n),x(:,2:n));
                D=Matrix(ind(1),:);
                %%
                class(ind(1))=no;
                ind(1)=[];
                i1=find(D<=Eps);              
                     
                if length(i1)>=k; %%%%
                    type(ob(1))=1;
                    for i=1:length(i1)
                        if touched(i1(i))==0
                            touched(i1(i))=1;
                            ind=[ind i1(i)];
                            class(i1(i))=no;
                        end
                    end
                else
                    type(ob(1))=0;
                end        %%%%%%%%%%%%%%%
                
            end
            no=no+1;
        end
    end
end

i1=find(class==0);
class(i1)=-1;
type(i1)=-1;
end
