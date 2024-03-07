function MIs  = FastPairMI(data,h)
% data : the input data, rows correspond to genes
%        columns correspond to arrays (samples)  
% h    : the std of the Gaussian kernel for density estimation 



MIs = zeros(size(data,1));
h_square = h^2;
L = size(data,2);
for k=1:L
    % tmp = data - repmat(data(:,k),1,L); same as the line below
    tmp = data - data(:,k);
    tmp = exp(-(tmp.^2)/(2*h_square));
    tmp1 = sum(tmp,2);

    % tmp2 = tmp*tmp';
    % for j=1:size(tmp2,1)
    %     tmp2(j,:) = tmp2(j,:)./tmp1(j);
    %     tmp2(:,j) = tmp2(:,j)./tmp1(j);
    % end
    % MIs = MIs + log(tmp2);
    %clear tmp2

% %   The following commented line does the same job as lines 16~22
     MIs = MIs + log((tmp*tmp')./(tmp1*tmp1'));
end
MIs = MIs/L + log(L);

end