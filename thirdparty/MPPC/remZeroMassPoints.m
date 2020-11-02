function [y,z,b,I,y_mass,cut_indices] = remZeroMassPoints(y,z,b,mass,I,cut_indices)
%remZeroPointsMass removes points on y with zero projected mass 

[m,~] = size(y);
y_mass = zeros(m,1); %the effective mass at each point on y, w.r.t. x_rad

for k=0:m-1
    
    I_k = I==m-k;
    y_mass(m-k) = sum(mass(I_k));
    
    if y_mass(m-k) == 0 %remove points with no mass
        [y,z,b] = removePoints(y,z,b,m-k);
        cut_indices = cut_indices - double(cut_indices>=m-k);
        y_mass(m-k) = [];
        I = I - (I>m-k);
    end
end

[m,~] = size(y);
cut_indices(cut_indices<1)=[];
cut_indices(cut_indices>m-1)=[];

end

