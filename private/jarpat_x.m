function [Fpn, conn]=jarpat_x(X,k,j,tr1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jarvis-Patrick clustering 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% X: data set
% k: numb. of neghbours taken into account
% j: min similarity considered as member of the same cluster
% tr1: min similarity considered as a neighbour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fpn=X;
            %Parameters
n=k;        %numb. of neghbours taken into account
tr1=tr1;	%min similarity considered as a neighbour
tr2=j/10;	%min similarity considered as member of the same cluster
            %Similarity matrix
% S=[];
% O=[];
% for i=1:size(Fpn,1)
%     dumm=[];
%     for j=1:size(Fpn,1)
%        dumm(j)=sqrt(sum((Fpn(i,:)-Fpn(j,:)).^2)); %calculating distances
%     end   
%     [dumm,I]=sort(dumm); 
%     S=[S;dumm(1:n)];
%     O=[O;I(1:n)];
% end   
% S=1-S/max(max(S));

[O,S]=knnsearch(X',X','K',k);
S=1-S/max(S(:));

% [O2,S2]=knnsearch(X',X','K',5);
% isequal(O,O2)
% S2=1-S2/max(S2(:));

%Clustering
conn=zeros(size(S,1),1); 	%Membership vector
cent=zeros(size(S,1),1);	%Indicator vector
j=0;
while max(conn==0)
   j=j+1;
   %disp(['Cal. of  the clusters ',num2str(j)])
   dumm=find(conn==0);	 
   start=dumm(1);       % The first in the new cluster
   cent(start)=-1;      % not tested yet
   conn(start)=j;       % this was the center
   
   while max(cent==-1) 
       dumm=find(cent==-1);	 
       start=dumm(1);   % The next untested center
       cent(start)=1;   % this was the center 
          
	   str1=O(start,S(start,:)>tr1);	 % The neighbours of this data
        B1=length(str1);
        TS=[];
        
   	% What is the sim. of this data in terms of neighbours  
   	for c=1:length(str1)
         str2=O(str1(c),S(str1(c),:)>tr1);
         
         B2=length(str2);
         dumm=sort([str1 str2]);
         BC=sum(dumm(1:end-1)==dumm(2:end));    % Numb. of the same neighb.
	     TS(c)=BC/(B1+B2-BC);                   % Tanimoto similarity
    end
        newconn=str1(TS>tr2);       % the new members of the cluster 
        conn(newconn)=j;
        cent(newconn(cent(newconn)==0))=-1; % these were not assigned (NEW)
          cent(start)=1;
          conn(start)=j;    % this was the center
   end
end

end
%text(Fpn(:,1),Fpn(:,2),num2str(conn))

