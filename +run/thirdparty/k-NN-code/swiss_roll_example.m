% Swiss roll example
%
%
% Written by Jose Costa (jcosta@umich.edu), Alfred Hero (hero@eecs.umich.edu), 2004
%
%
% BEGIN Copyright notice
% 
% Matab scripts for intrinsic dimension and entropy estimation using k-nearest
% neighbor graphs. The details of the algorithms can be found in:
% 
%   J. A. Costa and A. O Hero, "Entropic Graphs for Manifold Learning",
%   Proc. of IEEE Asilomar Conf. on Signals, Systems, and Computers,
%   Pacific Groove, CA, November, 2003.
%
%   J. A. Costa and A. O. Hero, "Geodesic Entropic Graphs for Dimension and
%   Entropy Estimation in Manifold Learning", 
%   to appear in IEEE Trans. on Signal Processing, Aug., 2004. 
%
% Published reports of research using the code provided here (or a modified version)
% should cite the two articles referenced above.
%
% Comments and questions are welcome. We would also appreciate hearing about 
% how you used this code, improvements made to it, etc. You are free to modify the
% code, as long as you reference the original contributors.
%
% END Copyright notice


n=400; % number of sample points

d=2;  % Intrinsic dimension

samp=(n-10):(n-1);% sample points to build bootstrap estimate of mean length function

kneighbors=4; % number of neighbors in kNN

gamma=1;




% generate swiss roll data
tt = (3*pi/2)*(1+2*rand(1,n));  height = 21*rand(1,n);

X = [tt.*cos(tt); height; tt.*sin(tt)];     % Matrix of points


%------------------------------------------------------
% k-NN graph examples ---------------------------------
%------------------------------------------------------

[dknn_1,Hknn_1,Lavgknn_1,Lstdknn_1]=knn_graph_estim_1(X',kneighbors,gamma,1,10,samp);       % (M,N)=(1,10)

display(['Intrinsic dimension estimate = ',num2str(round(dknn_1))])                         % note that dknn_1 is rounded to the nearest integer
display(['Intrinsic entropy estimate = ',num2str(Hknn_1),' bits'])



[dknn_2,Hknn_2,Lavgknn_2,Lstdknn_2]=knn_graph_estim_1(X',kneighbors,gamma,10,1,samp);       % (M,N)=(10,1)

display(['Intrinsic dimension estimate = ',num2str(round(dknn_2))])
display(['Intrinsic entropy estimate = ',num2str(Hknn_2),' bits'])




Dist = L2_distance(X,X,1);      % matrix of Euclidean distances between points
% Dumb implementation of k-NN algorithm. Much slower in low (extrinsic) dimensional spaces (like this example)
% but faster for high (extrinsic) dimensional spaces (like images)
[dknn_3,Hknn_3,Lavgknn_3,Lstdknn_3]=knn_graph_estim_2(Dist,kneighbors,gamma,1,10,samp);       % (M,N)=(1,10)

display(['Intrinsic dimension estimate = ',num2str(round(dknn_3))])
display(['Intrinsic entropy estimate = ',num2str(Hknn_3),' bits'])



%%

figure
colormap jet; %set(gcf,'Position',[200,400,620,200]);
subplot(1,2,1);
scatter3(X(1,:),X(2,:),X(3,:),12,'o', 'filled');
axis equal
view(-19.5,12)
title('Original Data');


subplot(1,2,2);
%cla;
plot_knn(X,kneighbors)
view(-19.5,12)
title('k-NN graph');
axis equal

%axis equal;