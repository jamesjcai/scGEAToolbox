function [ disc_energy,cont_energy ] = calculateEnergy( y,x,mass,lambda1,lambda2,I,cut_indices )
%CALCULATEENERGY Calculates both the discrete and continuous MPPC energies

if isempty(I)
    I = dsearchn(y,x)';
end
I_partial = computeProjs(x,y,I,cut_indices);

m = size(y(:,1));
normv = sqrt(sum((y(2:m,:)-y(1:m-1,:)).^2,2));
normv(cut_indices) = 0;
y_leng = sum(normv);
disc_energy = lambda1*(y_leng+lambda2*length(cut_indices)) + mass*sum((x-y(I,:)).^2,2);
cont_energy = lambda1*(y_leng+lambda2*length(cut_indices)) + mass*sum((x-getYVals(y,I_partial)).^2,2);

end


function yvals = getYVals(y,partial_ind)
[~,d] = size(y);
partials = mod(partial_ind,1)';
yvals = repmat((1-partials),1,d).*y(floor(partial_ind),:) + repmat(partials,1,d).*y(ceil(partial_ind),:);

end