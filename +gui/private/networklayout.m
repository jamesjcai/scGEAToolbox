function [x, y] = networklayout( W, nodeSizes, xratio, x, y)
    if(nargin < 3); xratio = 1; end
    xcenter = xratio/2;

    nodeSizes = nodeSizes / 100;
    nNode = size(W, 1);

    if(nargin < 4); x = rand(1, nNode); end
    if(nargin < 5); y = rand(1, nNode); end
    if(isempty(x)); x = rand(1, nNode); end
    if(isempty(y)); y = rand(1, nNode); end
    if(iscolumn(x)); x = x'; end
    if(iscolumn(y)); y = y'; end
    
    Dt = nodeSizes + nodeSizes';
    Dt = Dt - diag(diag(Dt));

    alpha = 0.1/nNode;

    w = full(W);
    numStep = 1000;
    stepSizes = logspace(-1, -3, numStep);
%     stepSizes = [stepSizes, 1e-3*ones(1,6000)];
%     numStep = length(stepSizes);
    for i = 1:numStep
        stepSize = stepSizes(i);
        x(x>=(1.02*xratio)) = xratio;
        y(y>=1.02) = 1.0;
        x(x<=(-0.02*xratio)) = -0.0;
        y(y<=-0.02) = -0.0;
        D = squareform(pdist([x/xratio;y]'));
        xcomp = xratio*(x - x')./D;
        ycomp = (y - y')./D;
        xcomp(isnan(xcomp)) = 0;
        ycomp(isnan(ycomp)) = 0;

%         overlap_coeff = 1 .* min(i, 50) / 50;
        overlap_coeff = 1;
        
        Dpull = D;
        Dpull = Dpull - Dt*1;
        overlapping_indices = D <= (Dt*overlap_coeff);
        Dpull(Dpull >= 0.6) = 0.6;
        Dpull(overlapping_indices) = 0;
        
        Dpush = alpha * (1./(Dpull + 0.2))/2;
        Dpush(overlapping_indices) = 0;

        Doverlap = zeros(size(D));
        Doverlap(overlapping_indices) = 1.2;

        F = stepSize * (Doverlap + Dpush .* ~w - Dpull .* w);
        Fx = F .* xcomp;
        Fy = F .* ycomp;
        Fx = Fx - diag(diag(Fx));
        Fy = Fy - diag(diag(Fy));    
        Fcx = stepSize * (xcenter - x).*(xcenter - x).*(xcenter - x) * 0.7 / (xratio^3);
        Fcy = stepSize * (0.5 - y).*(0.5 - y).*(0.5 - y) * 0.7;
        a = sum(Fx, 1) + Fcx;
        b = sum(Fy, 1) + Fcy;
        x = x + a;
        y = y + b;
    end
    x = x';
    y = y';
end

