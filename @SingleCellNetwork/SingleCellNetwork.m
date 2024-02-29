classdef SingleCellNetwork
    properties
        G % digraph
    end

    properties (Dependent)
        NumGenes
        A
        g
    end

    methods

        s = katz_score(obj, A, b)
        s = katz_centrality(obj, a, b)
        obj = filtadjc(obj, q)
        T = virtualknockout(obj, gid)
        [scnout] = subnetwork(obj, gid)
        [scnout] = subnetwork_neighbors(obj, targetg, siz)

        function obj = SingleCellNetwork(A, g)
            if nargin < 1, A = []; end
            if nargin < 2, g = []; end
            if ~isempty(A) && isempty(g)
                g = textscan(sprintf('g%d,', 1:size(A, 1)), '%s', 'delimiter', ',');
                g = string(g{1});
                %g=string(transpose(1:size(A,1)));
                %for k=1:size(A,1)
                %    g(k)=sprintf("g%s",g(k));
                %end
            end
            if ~isempty(g)
                g = matlab.lang.makeUniqueStrings(g);
            end
            obj.G = digraph(A, g, 'omitselfloops');
        end

        function A = get.A(obj)
            A = adjacency(obj.G, 'weighted');
        end

        function g = get.g(obj)
            g = string(obj.G.Nodes.Variables);
        end

        function m = get.NumGenes(obj)
            m = size(obj.G.Nodes, 1);
        end

        function obj = set.NumGenes(obj, ~)
            fprintf('%s%d\n', 'NumGenes is: ', obj.NumGenes)
            error('You cannot set NumGenes property');
        end

        function grnview(obj)
            sc_grnview(obj.A, obj.g);
        end

        function p = plot(obj, allgenes)
            if nargin < 2 && obj.NumGenes <= 500
                allgenes = true;
            elseif nargin < 2 && obj.NumGenes > 500
                allgenes = false;
                disp('Cannot display >500 nodes.');
                disp('A random subgraph is shown.')
                disp('To force to display all nodes, use plot(SCN,true);');
            end
            if allgenes
                p = plot(obj.G);
            else
                rng("shuffle");
                rid = randperm(obj.NumGenes);
                xg = subgraph(obj.G, rid(1:100));
                [bin, binsize] = conncomp(xg, 'Type', 'weak');
                idx = binsize(bin) == max(binsize);
                SG = subgraph(xg, idx);
                p = plot(SG);
                title('Random Subnetwork')
            end
        end

        function p = plotweighted(obj)
            % if nargin < 2 && obj.NumGenes <= 500
            %     allgenes = true;
            % elseif nargin < 2 && obj.NumGenes > 500
            %     allgenes = false;
            %     disp('Cannot display >500 nodes.');
            %     disp('A random subgraph is shown.')
            %     disp('To force to display all nodes, use plot(SCN,true);');
            % end
            % if allgenes
            %     p = plot(obj.G);
            %     obj.G.Edges.LWidths = 7 * obj.G.Edges.Weight / max(abs(obj.G.Edges.Weight));
            %     p.LineWidth = abs(obj.G.Edges.LWidths);
            % else
            %     rng("shuffle");
            %     rid = randperm(obj.NumGenes);
            %     xg = subgraph(obj.G, rid(1:100));
            %     [bin, binsize] = conncomp(xg, 'Type', 'weak');
            %     idx = binsize(bin) == max(binsize);
            %     SG = subgraph(xg, idx);
            %     %p=plotweighted(SG);
            %     LWidths = 7 * SG.Edges.Weight / max(SG.Edges.Weight);
            %     p = plot(SG, 'LineWidth', abs(LWidths));
            %     title('Random Subnetwork')
            % end
            sc_grnview(obj.A,obj.g);
        end


        function obj = makesparse(obj, q)
            if nargin < 2, q = 0.95; end
            a = max(abs(obj.A(:)));
            if a == 0, a = 1; end
            obj.A = obj.A ./ a;
            obj.A = obj.A .* (abs(obj.A) > quantile(abs(obj.A(:)), q));
            if ~issparse(obj.A)
                obj.A = sparse(obj.A);
            end
        end

    end
end
