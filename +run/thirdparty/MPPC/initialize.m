function [ y,cut_indices,energy ] = initialize(x,mass,lambda1,lambda2 )
%INITIALIZE initializes y to a set of singletons
%   The initialization consists of a two-step line search for the number of
%   singletons that minimizes the (MPPC) functional.

n = length(x(:,1));
energy_prev = inf;
k = 1;
[I_new,y_new] = kmeans(x,k); cut_indices = (1:k-1)';
energy = calculateEnergy(y_new,x,mass,lambda1,lambda2,I_new',cut_indices);
warning('off','stats:kmeans:FailedToConverge')
while energy < energy_prev && k<n
    y = y_new;
    energy_prev = energy;
    k = min(2*k,n);
    [I_new,y_new] = kmeans(x,k,'MaxIter',100); cut_indices = (1:k-1)';
    energy = calculateEnergy(y_new,x,mass,lambda1,lambda2,I_new',cut_indices);
end
if energy<energy_prev
    y = y_new;
end
[I,y] = kmeans(x,[],'Start',y); cut_indices = (1:length(y(:,1))-1)';
energy_prev = calculateEnergy(y,x,mass,lambda1,lambda2,I',cut_indices);

while k>1 && length(y(:,1))<n
    k = min(floor(k/2),n-length(y(:,1)));
    y_new = addSingletons(x,mass,I,y,k);
    [I_new,y_new] = kmeans(x,[],'Start',y_new,'MaxIter',100); cut_indices = (1:length(y_new(:,1))-1)';
    energy = calculateEnergy(y_new,x,mass,lambda1,lambda2,I_new',cut_indices);
    if energy < energy_prev
        y = y_new; I = I_new;
        energy_prev = energy;
    end
end
[I,y] = kmeans(x,[],'Start',y); cut_indices = (1:length(y(:,1))-1)';
z = y(2:end,:)-y(1:end-1,:) ; b = zeros(size(z));
[y,~,~,cut_indices] = checkSingletons(x,mass,I,y,z,b,cut_indices,lambda1,lambda2);
energy = calculateEnergy(y,x,mass,lambda1,lambda2,[],cut_indices);

end


function [ y ] = addSingletons( x,mass,I,y,k)
% This function adds k singletons to y near the existing singletons where 
% the data fidelity is the highest

[m,d] = size(y);
fid = zeros(m,1);

for i=1:m
    x_ind = I == i;
    n1 = sum(x_ind);
    fid(i) = mass(x_ind)*sum((repmat(y(i,:),n1,1)-x(x_ind,:)).^2,2);
end

[~,sorti] = sort(fid,'descend');
new_points = zeros(k,d);

for i=1:k
    x1 = x(I==sorti(i),:); 
    % location for new singleton is between existing singleton and a data
    % point that projects to it
    new_points(i,:) = .5*(x1(1,:) + y(sorti(i),:)); 
end
y = [y;new_points];

end

