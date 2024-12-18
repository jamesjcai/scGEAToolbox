function i_networkvis(G, xy, curved, fontsize)

if nargin < 3, curved = true; end
if nargin < 4, fontsize = 15; end

textOpts.FontSize = fontsize;
textOpts.HorizontalAlignment = 'center';
textOpts.VerticalAlignment = 'middle';
textOpts.FontWeight = 'normal';

hFig = figure("Visible","off", 'SizeChangedFcn', @in_addtext);

%if ~curved
%    gplot(G.adjacency, xy, '-k');
%else
    curve_gplot(G.adjacency, xy, curved);
%end
    
    hold on
    scatter(xy(:,1), xy(:,2), 300, ...
        'MarkerEdgeColor','k', ...
        'MarkerFaceColor',[.8 .8 .8]);

    % for k=1:length(G.Nodes.Name)
    %     [wx, hx] = measureText(G.Nodes.Name{k}, textOpts);
    %     text(xy(k,1)-wx/2, xy(k,2), ...
    %         G.Nodes.Name{k},'FontSize',15,...
    %         'BackgroundColor','w', ...
    %         'FontWeight','normal', ...
    %         'HorizontalAlignment','center', ...
    %         'VerticalAlignment','middle','Margin',0.2);
    % end
    in_addtext([],[]);
    hFig.Visible='on';


    function in_addtext(~, ~)
        textHandles = findall(gcf, 'Type', 'text'); % Find all text objects in the current figure
        delete(textHandles); % Delete all found text objects        
        for k=1:length(G.Nodes.Name)
            [wx, hx] = measureText(G.Nodes.Name{k}, textOpts);
            text(xy(k,1)-floor(wx/2), xy(k,2), ...
                G.Nodes.Name{k},'FontSize',textOpts.FontSize,...
                'BackgroundColor','w', ...
                'FontWeight','normal', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','Margin',0.2);
        end
    end

end


 function [width, height] = measureText(txt, textOpts, axis)
    if(nargin < 3)
       axis = gca(); 
    end
    if nargin < 2
        textOpts = struct();
        textOpts.HorizontalAlignment = 'center';
        textOpts.VerticalAlignment = 'middle';
        textOpts.FontSize = 20;
        textOpts.FontWeight = 'normal';
    end
    hTest = text(axis, 0, 0, txt, textOpts);
    textExt = get(hTest, 'Extent');
    delete(hTest);
    height = textExt(4)/3;    %Height
    width = textExt(3)/3;     %Width
 end


 function curve_gplot(sbeG, xy, curved)
 if nargin<3, curved = false; end

    [i, j] = find(sbeG);
    [~, p] = sort(max(i, j));
    i = i(p);
    j = j(p);
    
    X = [xy(i, 1), xy(j, 1)]';
    Y = [xy(i, 2), xy(j, 2)]';
    
    if isfloat(xy)
        X = [X; NaN(size(i))'];
        Y = [Y; NaN(size(i))'];
    end
    X = X(:);
    Y = Y(:);
    
    
    M = zeros(1, 4);
    hold on
    for k = 1:length(X) - 1
        if ~(isnan(X(k)) || isnan(X(k+1)))
            if curved
                p1 = [X(k), Y(k)];
                p2 = [X(k+1), Y(k+1)];
        
                if ~ismember([p1, p2], M, 'rows')
                    [a, b] = quadraticcurveto(p1, p2);
                    plot(a, b, '-', ...
                        'color', [.6, .6, .6], 'linewidth', 1);
                    
                    M = [M; [p1, p2]];
                    M = [M; [p2, p1]];
                end
            else
                
                plot([X(k),X(k+1)],[Y(k),Y(k+1)],'r:');
            end
        end
    end
 end


 function [x, y] = quadraticcurveto(currentp, curvetop, controllp)

    %The quadraticCurveTo method creates a line from the path's current point to the specified point, via a controlpoint.
    %The quadraticCurveTo method takes four parameters, the x and y coordinates for the controlpoint, and the x and y coordinates for the line's destination.

    %Current point moveTo(20,20)
    %CurveTo point quadraticCurveTo(20,100,200,20)
    %Controllpoint quadraticCurveTo(20,100,200,20)

    o = curvetop(1);
    i = curvetop(2);
    d = currentp(1);
    f = currentp(2);

    x2 = curvetop(1);
    %y2=curvetop(2);
    x1 = currentp(1);
    y1 = currentp(2);

    if nargin < 3
        %x3=min([x1,x2])+abs(x1-x2)/2;
        %y3=y1;    controllp=[x3 y3];
        controllp = [(d + o) / 2 + (i - f) / 4, (f + i) / 2 + (d - o) / 4];
    end
    P = [[currentp'; 0], [controllp'; 0], [curvetop'; 0]];

    n = 3;
    count = 1;
    div = 50; %number of segments of the curve (Increase this value to obtain a
    %smoother curve
    for u = 0:(1 / div):1
        sum = [0, 0, 0]';
        for i = 1:n
            B = nchoosek(n, i-1) * (u^(i - 1)) * ((1 - u)^(n - i + 1)); %B is the Bernstein polynomial value
            sum = sum + B * P(:, i);
        end
        B = nchoosek(n, n) * (u^(n));
        sum = sum + B * P(:, n);
        A(:, count) = sum; %the matrix containing the points of curve as column vectors.
        count = count + 1; % count is the index of the points on the curve.
    end

    x = A(1, :);
    y = A(2, :);
 end