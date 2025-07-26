function i_plot_pseudotimeseries(X, genelist, t, genes)
    % Plot pseudotime series for selected genes
    
    % Default genes to plot if not provided
    if nargin < 4
        genes = string(['AICDA', 'BACH2', 'BCL6', 'IRF4', 'PAX5', 'PRDM1', 'REL', 'RELA']);
    end
    
    % Sort t (pseudotime) and re-order X accordingly
    [t, i] = sort(t);
    
    % Get indices of the genes of interest
    n = length(genes);
    idx = [];
    glist = [];
    for k = 1:n
        ixx = find(genelist == genes(k));
        if ~isempty(ixx)
            idx = [idx; ixx];
            glist = [glist; genes(k)];
        end
    end
    
    % Subset the expression data
    x = X(idx, i);
    
    % Define colors for plotting
    co = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; ...
          0.6350, 0.0780, 0.1840; 0, 0, 1];
    co = repmat(co, ceil(length(idx)/size(co,1)), 1); % Extend color array if needed
    set(groot, 'defaultAxesColorOrder', co(1:length(idx), :));
    
    % Plot the pseudotime series
    % figure;
    hold on;
    msz = 1;
    for k = 1:size(x, 1)
        plot(t, x(k, :), '.', 'MarkerSize', msz);
    end
    xlabel('Pseudotime');
    ylabel('Expression');
    
    % Add locfit and predict paths
    pw1 = fileparts(mfilename('fullpath'));
    locfit_paths = {fullfile(pw1, '..', 'external', 'fun_locfit')};
    
    for p = locfit_paths
        if ~(ismcc || isdeployed)
            addpath(p{1});
        end
    end
    
    % Fit smooth curves using locfit
    if size(t, 2) ~= 1
        t = t'; % Ensure t is a column vector
    end
    
    Pk = [];
    for k = 1:size(x, 1)
        fitm1 = locfit(t, x(k, :)');  % Fit pseudotime series
        Pk = [Pk, plot(t, predict(fitm1, t), '-', 'LineWidth', 3)];
    end
    
    % Add legend and final touches
    box on;
    legend(Pk, glist, 'Location', 'eastoutside');
    xlim([0, max(t)]);
    
end
