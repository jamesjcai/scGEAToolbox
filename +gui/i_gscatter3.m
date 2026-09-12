function [h] = i_gscatter3(s, c, methodid, targetc, hAx)

if nargin < 5, hAx = []; end
if nargin < 4, targetc = 1; end
if nargin < 3, methodid = 1; end
if nargin < 2, c = ones(size(s, 1), 1); end

x = s(:, 1);
y = s(:, 2);

% STRING() throws on a sparse array ("Conversion to string from sparse
% single is not possible"), and SCE.X is single sparse from R2025a on, so
% anything derived from it by summing - library size, detected genes, the
% mt ratio - arrives here sparse. C is also handed to SCATTER as CData
% below, which rejects sparse too. Normalise once here rather than
% trusting every caller: a grouping variable is never usefully sparse.
% Non-numeric C (cellstr, string, categorical) is left alone for
% FINDGROUPS.
%
% FINDGROUPS is applied to the NUMBERS, not to STRING(c). Routing numeric
% C through STRING() ordered the groups lexicographically, so from ten
% groups on the colours were scrambled against the values: with clusters
% 1..12, cluster 10 took colour index 2 and cluster 2 took index 5. That
% also silently broke the caller's colorbar, because scgeatoolApp builds
% its tick labels from its own FINDGROUPS and labels colour index k with
% cL(k) - correct only while this function preserves that order.
if isnumeric(c) || islogical(c)
    c = full(double(c(:)));
    c = findgroups(c);
    % FINDGROUPS returns NaN for NaN input. Under the old STRING() path
    % those became an ordinary "NaN" group with a colour of its own;
    % leaving them NaN would instead drop the points from the plot without
    % saying so, so keep them as one group, placed last.
    isn = isnan(c);
    if any(isn)
        c(isn) = max([0; c(~isn)]) + 1;
    end
else
    c = findgroups(string(c));
end
kc = numel(unique(c));

if size(s, 2) >= 3, z = s(:, 3); end

switch methodid
    case 1
        if size(s, 2) == 2
            if isempty(hAx)
                h = scatter(x, y, 10, c);
            else
                h = scatter(hAx, x, y, 10, c);
            end
        elseif size(s, 2) >= 3

            if isempty(hAx)
                h = scatter3(x, y, z, 10, c);
            else
                % hAx
                % assignin("base", "c", c);
                % assignin("base", "x", x);
                % assignin("base", "y", y);
                % assignin("base", "z", z);
                assert(numel(c)==numel(x));
                h = scatter3(hAx, x, y, z, 10, c);
            end
        end
    case 2
        if size(s, 2) == 2
            h = gscatter(x, y, c, [], [], 5, 'off');
        elseif size(s, 2) >= 3
            h = gscatter3b(x, y, z, c, [], [], 5, 'off');
        end
        box off
    case 3
        h = gui.i_gscatter3(s, c, 1);
        h.MarkerEdgeAlpha = 0;
        hold on
        idx = c == targetc;
        h = gui.i_gscatter3(s(idx, :), c(idx), 1);
        hold off
end


if isempty(hAx)
        grid on
        colormap(pkg.i_mycolorlines(kc));
    else
        grid(hAx,'on');
        colormap(hAx,pkg.i_mycolorlines(kc));
    end
end
