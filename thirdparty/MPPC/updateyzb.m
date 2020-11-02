function [ y_new,z_new,b_new,I,energy,cut_indices ] = updateyzb(y,z,b,x,mass,lambda,rho,cut_indices,I)
%UPDATEYZB
% This function performs one step of the ADMM/split bregman algorithm, updating
% each of y,z,b once. It also takes into account seperate components of y
% given by cut_indices. If multiple components exist, the ADMM step is run 
% once on each component with its corresonding data and mass.
% If a particular component has no mass (no corresponding data), then it is
% deleted entirely from y, and z,b are adjusted accordingly.

% INPUTS:
% y-current points for curve
% z-variable for split bregman
% b-variable for split bregman
% x-the data points
% mass-the mass of the data points
% lambda-the coefficient of the length penalty
% rho-coefficient of the penalty for the discrepancy between phi(y) and z
% cut_indices-the last y indices of each component
% I - if it is nonempty, this one will be used (new projections won't be
% computed)
%
% OUTPUTS:
% I-the assignments of data to closest points in y
% energy-the energy functional evaluated at y_new, not including the
% lambda2 term
% the rest are updated versions of the corresponding inputs


[m,d] = size(y);
z0 = y(2:end,:)-y(1:end-1,:);

%Step 1: update z's
if isempty(z) %for the first iteration
    z = z0;
else
    prez = z0+b;
    normsprez = sqrt(sum(prez.^2,2)); % computes the euclidean norm of each row vector
    normalprez = bsxfun(@rdivide,prez, normsprez);
    z = bsxfun(@times,normalprez, max(normsprez-lambda/rho,0));
end

%Step 2: update b's
if isempty(b) %for the first iteration
    b = zeros(m-1,d);
else
    b = b + z0 - z;
end

%Step 3: update y's
if isempty(I)
    [I,~] = dsearchn(y,x);
    I = I';
end

[y,z_new,b_new,I,y_mass,cut_indices] = remZeroMassPoints(y,z,b,mass,I,cut_indices);
[m,~] = size(y);

if isempty(cut_indices)
    y_new = updateY(y,z_new,b_new,x,I,mass,y_mass,rho);
    energy = mass*sum((x-y_new(I,:)).^2,2);
    if m>1
        energy = energy + lambda*sum(sqrt(sum((y_new(2:end,:)-y_new(1:end-1,:)).^2,2)));
    end
    
else
    y_new = zeros(size(y));
    cut_indices = union(cut_indices(1),cut_indices); %to avoid repeated values
    cut_partition = [0;cut_indices;m];
    
    energy1 = 0;
    
    for j=2:length(cut_indices)+2   %perform updates on each component. This is parallelizable
        
        comp_y_indices = cut_partition(j-1)+1:cut_partition(j);
        comp_y = y(comp_y_indices,:);
        comp_ymass = y_mass(comp_y_indices);
        
        comp_xindices = find((I>cut_partition(j-1) & I<=cut_partition(j)));
        comp_x = x(comp_xindices,:);
        comp_mass = mass(comp_xindices);
        
        comp_I = I(comp_xindices)-cut_partition(j-1)*ones(1,length(comp_xindices));
        comp_z_indices = cut_partition(j-1)+1:cut_partition(j)-1;
        comp_z = z_new(comp_z_indices,:);
        comp_b = b_new(comp_z_indices,:);
        
        comp_y = updateY(comp_y,comp_z,comp_b,comp_x,...
            comp_I,comp_mass,comp_ymass,rho);
        
        y_new(comp_y_indices,:) = comp_y;
        
        if length(comp_y(:,1)) > 1
            energy1 = energy1 + lambda*sum(sqrt(sum((comp_y(2:end,:)-comp_y(1:end-1,:)).^2,2)));
        end
        
    end
    energy = energy1 + mass*sum((x-y_new(I,:)).^2,2);
end
end




function y_new = updateY(y,z,b,x,I,mass,y_mass,rho)
%Helper function used for updating y on a single component.

[m,d] = size(y);

x_bar = zeros(m,d);
for i=1:m
    assert(y_mass(i)~=0);
    x_ind = I==i;
    x_bar(i,:) = mass(x_ind)*x(x_ind,:);
end

if m==1 % If singleton, then position at mean of data projecting to it
    y_new = bsxfun(@rdivide,x_bar,y_mass);
else
    D = 2*y_mass+2*rho*ones(m,1);
    D(1) = 2*y_mass(1) + rho;
    D(m) = 2*y_mass(m) + rho;
    zb1 = [zeros(1,d);z-b];
    zb2 = [z-b;zeros(1,d)];
    zb = zb1 - zb2;
    RHS = 2*x_bar + rho*zb;
    y_new = tridiag((-rho)*ones(m-1,1),D,(-rho)*ones(m-1,1),RHS);
end

end


