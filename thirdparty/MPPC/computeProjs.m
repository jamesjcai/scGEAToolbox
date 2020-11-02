function [ I_new ] = computeProjs( x,y,I,cut_indices )
%computeProjs computes non-discrete projections to the polygonal curve, 
%given the closest point in y.
% This function computes the projections onto edges of y, and is used in
% computing the continuous MPPC energy. Note that projections computed here
% are approximate, given projections onto the vertices of y. That is, this 
% function will only project data to the edge (y_i, y_{i+1}) if the data
% either projects to y_i or y_{i+1}.
% The argument 'I' may be passed as empty, in which case it will be
% computed here.

[m,d] = size(y);
[n,~] = size(x);

if isempty(I)
    I = dsearchn(y,x)';
end
I_new = I;

v = y(2:m,:)-y(1:m-1,:);
normv = sqrt(sum(v.^2,2));

start_indices = [1;cut_indices+1];
end_indices = [cut_indices;m];

singleton_ind = start_indices==end_indices;
singletons = start_indices(singleton_ind);

start_indices = start_indices(~singleton_ind);
end_indices = end_indices(~singleton_ind);

% find x indices that are closest to first vertex of component
x_aff1 = ismember(I,start_indices);
x_aff1i = find(x_aff1);
if ~isempty(x_aff1i)
    % indices that project onto first edge of component
    x_aff = x_aff1i((dot((y(I(x_aff1i),:)-x(x_aff1i,:))',v(I(x_aff1i),:)') < 0));
    %update I-value. It will be non-integer (e.g. 1.5 would mean it projects to
    %the midpoint between y(1) and y(2))
    if ~isempty(x_aff)
        I_new(x_aff) = I(x_aff) + .5 * ( 1 + ( sum((y(I(x_aff),:)-x(x_aff,:)).^2,2) - ...
            sum((y(I(x_aff)+1,:)-x(x_aff,:)).^2,2) ) ./ (sum((y(I(x_aff),:)-y(I(x_aff)+1,:)).^2,2)) )';
    end
    
end
% find x indices that are closest to last vertex of component
x_aff2 = ismember(I,end_indices);
x_aff2i = find(x_aff2);
if ~isempty(x_aff2i)
    % indices that project onto last edge of component
    x_aff = x_aff2i((dot((y(I(x_aff2i),:)-x(x_aff2i,:))',v(I(x_aff2i)-1,:)') > 0));
    %update I-value.
    if ~isempty(x_aff)
        I_new(x_aff) = I(x_aff)-1 + .5 * ( 1 + ( sum((y(I(x_aff)-1,:)-x(x_aff,:)).^2,2)- ...
            sum((y(I(x_aff),:)-x(x_aff,:)).^2,2) ) ./ (sum((y(I(x_aff)-1,:)-y(I(x_aff),:)).^2,2)) )';
    end
end

% find x indices that are closest to interior vertex of component
x_aff = ~x_aff1 & ~x_aff2 & ~ismember(I,singletons);
x_affi = find(x_aff);
if ~isempty(x_affi)
    % indices that project onto interior edge of component
    unit_disp = (y(I(x_aff),:)-x(x_aff,:))./repmat(sqrt(sum((x(x_aff,:)-y(I(x_aff),:)).^2,2)),1,d);
    left_unit = v(I(x_aff)-1,:)./repmat(normv(I(x_aff)-1,:),1,d);
    right_unit = -v(I(x_aff),:)./repmat(normv(I(x_aff),:),1,d);
    left_dot = dot(left_unit',unit_disp');
    right_dot = dot(right_unit',unit_disp');
    % indices that project onto right edge
    x_aff_r = x_affi((right_dot > 0 & right_dot>left_dot));
    if ~isempty(x_aff_r)
        %update I value.
        I_new(x_aff_r) = I(x_aff_r) + .5 * ( 1 + ( sum((y(I(x_aff_r),:)-x(x_aff_r,:)).^2,2) - ...
            sum((y(I(x_aff_r)+1,:)-x(x_aff_r,:)).^2,2) ) ./ (sum((y(I(x_aff_r),:)-y(I(x_aff_r)+1,:)).^2,2)) )';
    end
    % indices that project onto left edge
    x_aff_l = x_affi((left_dot > 0 & left_dot>right_dot));
    if ~isempty(x_aff_l)
        %update I value.
        I_new(x_aff_l) = I(x_aff_l)-1 + .5 * ( 1 + ( sum((y(I(x_aff_l)-1,:)-x(x_aff_l,:)).^2,2)- ...
            sum((y(I(x_aff_l),:)-x(x_aff_l,:)).^2,2) ) ./ (sum((y(I(x_aff_l)-1,:)-y(I(x_aff_l),:)).^2,2)) )';
    end
    
end
I_new(I_new>m) = m; % account for numerical imprecision
I_new(I_new<1) = 1;


end

