%{
figure;
load ethanol;         % load the dataset.
   fit = locfit(E,NOx);   % local regression, with x,y vectors.
   
subplot(2,2,1)   
   lfplot(fit)           % plot the fitted curve.


X=E;
y=NOx;

   [smoothed_y, sortIdx] = loess_smoothing(X, y, 0.65);

   subplot(2,2,2)   

scatter(E, NOx);
hold on
plot(X(sortIdx), smoothed_y);


function [smoothed_y, sortIdx] = loess_smoothing(X, y, span)
    % loess_smoothing fits/smooths response data (y)
    % MATLAB documentation https://www.mathworks.com/help/curvefit/smooth.html
    %
    % INPUT:
    % X ============> Count matrix (column) for ig gene
    % y ============> Dependent function (e.g. spline speudotime)
    % span =========> Smoothing parameter, 0.3 preserves data noise, and 
    %                 0.75 makes it a curve
    %
    % OUTPUT:  
    % smoothed_y ===> Smoothed data corresponding to y
    % sortIdx ======> New sort for cells and response y
    %
    % USAGE:
    %
    % [y,idx] = ismember('splinefit_pseudotime', sce.list_cell_attributes(1:2:end));
    % assert(all(y))
    % t = sce.list_cell_attributes{idx+1};
    % ngene = size(X_qubo, 1);
    % ncell = size(X_qubo, 2);
    % y_qubo_fit = zeros(ncell, ngene);
    % t_qubo_sort = zeros(ncell, ngene);
    % % Loop over genes
    % for ig = 1:ngene
    %     [y_qubo_fit(1:ncell, ig), idx] = ...
    %                         loess_smoothing(t, X_qubo(ig,:)', 0.75);
    %     t_qubo_sort(1:ncell, ig) = t(idx);
    % end

    % Ensure inputs are column vectors
    X = X(:);
    y = y(:);

    % Sort the data by X
    [X_sorted, sortIdx] = sort(X);
    y_sorted = y(sortIdx);

    % method options: 
    % 'loess' Local regression - weighted linear least squares w 2nd deg pol.
    % 'lowess' Local regression - weighted linear least squares w 1st deg pol.
    % 'rloess' Robust ver of 'loess'
    % 'rlowess' Robust ver of 'lowess' but more expensive (not much improvement)
    method = 'loess'; %Very efficient and good prediction

    % Apply LOESS smoothing second degree polynomial
    smoothed_y = smooth(X_sorted, y_sorted, span, method);

    %smoothed_y = smooth(X_sorted, y_sorted, span, 'lowess');

    % %Plot the results
    % figure;
    % plot(X_sorted, y_sorted, 'bo', 'DisplayName', 'Original Data'); hold on;
    % plot(X_sorted, smoothed_y, 'r-', 'LineWidth', 2, 'DisplayName', 'LOESS Smoothed Data');
    % legend('show');
    % title(strcat(method, ' Smoothing'));
    % xlabel('X');
    % ylabel('y');
end
%}