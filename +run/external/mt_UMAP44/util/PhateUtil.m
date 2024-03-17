classdef PhateUtil < handle
%PHATEUTIL is a wrapper for the algorithm created at the Krishnaswamy lab.
%   The algorithm is described at
%   https://www.krishnaswamylab.org/projects/phate
%
%   This file is a wrapper that 
%   - Adds checking for PHATE inputs using MATLAB's addParameter idiom.  
%   - Provides meta info for our Args.m module, which will parse the
%     comments in this file and presents an "argument clinic" GUI.
% 
%   Most of the remaining comments are a direct copy from their phate.m
%   file. You can download it at https://github.com/KrishnaswamyLab/PHATE.
%
%   Y = phate(data) runs PHATE on data (rows: samples, columns: features)
%   with default parameter settings and returns a 2 dimensional embedding.
%
%   If data is sparse PCA without mean centering will be done to maintain
%   low memory footprint. If data is dense then normal PCA (with mean
%   centering) is done.
%
%
%   REQUIRED INPUT ARGUMENT
%   The argument csv_file_or_data is either 
%   A) a char array identifying a CSV text file containing the data to be reduced.
%   B) a numeric matrix to be reduced; a numeric matrix.
%
%   OPTIONAL NAME VALUE PAIR ARGUMENTS
%   The optional argument name/value pairs are:
% 
%   NAME            VALUE
% 
%   'ndim'          Number of (output) embedding dimensions. Common values
%                   are 2 or 3.
%                   Default is 2.
%
%   'k'             Number of nearest neighbors for bandwidth of adaptive
%                   alpha decaying kernel or, when a=[], number of nearest
%                   neighbors of the knn graph. For the unweighted kernel,
%                   we recommend k to be a bit larger, e.g. 10 or 15. For
%                   flow cytometry, we recommend 15.
%                   Default is 15.
%
%   'a'             Alpha of alpha decaying kernel. When a=[] or 0, knn 
%                   (unweighted) kernel is used. 
%                   Default is 40.
%
%   't'             Number of diffusion steps. 0 or [] automatically 
%                   picks the optimal t.
%                   Default is [] .
%
%   't_max'         Maximum t for finding optimal t. If t = [], optimal 
%                   t will be computed by computing Von Neumann Entropy 
%                   for each t <= t_max and then picking the kneepoint. 
%                   Default is 100.
%
%   'npca'          Number of PCA components for computing distances. 
%                   Default is 100.
%
%   'mds_method'    Method of multidimensional scaling. Choices are:
%                   - mmds - metric MDS 
%                   - cmds - classical MDS
%                   - nmmds - non-metric MDS
%                   Default is mmds
%
%   'distfun'       Distance function. 
%                   Default is euclidean.
%
%   'distfun_mds'   Distance function for MDS. 
%                   Default is euclidean.
%
%   'pot_method'    Method of computing the PHATE potential distance. 
%                   Choices are:
%                       - log: -log(P + eps). 
%                       - sqrt: sqrt(P). (not default but often produces 
%                           superior embeddings)
%                       - gamma: 2/(1-\gamma)*P^((1-\gamma)/2)
%                   Default is log.
%
%   'gamma'         Gamma value for gamma potential method. Value between
%                   -1 and 1. -1 is diffusion distance. 1 is log potential.
%                   0 is sqrt. Smaller gamma is a more locally sensitive
%                   embedding whereas larger gamma is a more globally
%                   sensitive embedding.
%                   Default is 0.5.
%
%   'pot_eps'       Epsilon value added to diffusion operator prior to
%                   computing potential. Only used when 'pot_method' is
%                   'log', i.e.: -log(P + pot_eps). 
%                   Default is 1e-7.
%
%   'n_landmarks'   Number of landmarks for fast and scalable PHATE. [] or
%                   n_landmarks = npoints does no landmarking, which is
%                   slower. More landmarks is more accurate, but comes at
%                   the cost of speed and memory.  For flow cytometry, we
%                   recommend 600.
%                   Default is 2000.
%
%   'nsvd'          Number of singular vectors for spectral clustering (for
%                   computing landmarks). 
%                   Default is 100.
%
%   'kernel'        User supplied kernel. If not given ([]) kernel is
%                   computed from the supplied data. Supplied kernel should
%                   be a square (samples by samples) symmetric affinity
%                   matrix. If kernel is supplied input data can be empty
%                   ([]).
%                   Defaults to [].

%   AUTHORSHIP
%   PHATE is the invention of Smita Krishnaswamy at Yale University.
%   The authors of this mere wrapper are
%   Primary Developer & Math Lead: Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
    
    properties(Constant)
        METRICS = {'euclidean';  'cityblock'; 'chebychev'; 'minkowski';...
            'mahalanobis'; 'seuclidean'; 'cosine'; 'correlation'; ...
            'spearman'; 'hamming'; 'jaccard'};
        METRIC_DICT=containers.Map({'euclidean', 'l2', 'manhattan', 'l1',...
            'taxicab', 'cityblock', 'seuclidean', 'standardised_euclidean',...
            'chebychev', 'linfinity', 'linfty', 'linf', 'minkowski',...
            'mahalanobis', 'cosine', 'correlation', 'hamming', 'jaccard',...
            'spearman'}, ...
            {'euclidean', 'euclidean', 'cityblock', 'cityblock',...
            'cityblock', 'cityblock', 'seuclidean', 'seuclidean', ...
            'chebychev','chebychev', 'chebychev', 'chebychev', 'minkowski', ...
            'mahalanobis','cosine', 'correlation', 'hamming', 'jaccard', ...
            'spearman'});
        POT_METHODS={'log', 'sqrt', 'gamma'};
        MDS_METHODS={'mmds', 'cmds', 'nmmds'};
        ARG_DICT=containers.Map({'ndim', 'n_components', 'k', 'n_neighbors',...
            'a', 't', 't_max', 'npca','mds_method','distfun','metric',...
            'distance', 'distfun_mds', 'pot_method', 'gamma', 'pot_eps',...
            'n_landmarks', 'nsvd', 'kernel'}, ...
            {'ndim', 'ndim', 'k', 'k',...
            'a', 't', 't_max', 'npca', 'mds_method','distfun', 'distfun'...
            'distfun', 'distfun_mds', 'pot_method', 'gamma', 'pot_eps',...
            'n_landmarks', 'nsvd', 'kernel'});
    end
    
    methods(Static)
        function args=Handle0(args)
            if ~isempty(args.t) && args.t==0
                args.t=[];%automatically find best
            end
            if ~isempty(args.a) && args.a==0
                args.a=[];%automatically find best
            end
        end

        function ok=IsInstalled
            try
                phate(randi(100, 10,5), 't', 2, ...
                    'progress_callback', 'none');
                ok=true;
                return;
            catch
                addpath ~/Downloads/PHATE-master/Matlab/
            end
            try
                phate(randi(100, 10,5), 't', 2);
                ok=true;
            catch ex
                ex.getReport
                ok=false;
                if askYesOrNo(Html.Sprintf(['Calling PHATE returned:' ...
                        '<br>"<i>%s</i>"<br><br><center>PHATE is ' ...
                        'provided by the Krishnaswamy Lab...<br>'...
                        '<b>Download PHATE from GitHub?</b></center><hr>'], ...
                        ex.message))
                    web(...
                        'https://github.com/KrishnaswamyLab/PHATE', ...
                        '-browser');
                    MatBasics.RunLater(@(h,e)explain, 1.5);
                end
            end

            function explain
                msg(Html.Sprintf( ...
                    ['<ul><li>Click the <font bgcolor="green" color="white">Code</font> button to download.' ...
                    '<li>Unzip the download.' ...
                    '<li>Call addpath on phate folder.<br>For example' ...
                    '<br>&nbsp;&nbsp;' ...
                    '<i>addpath ~/Downloads/PHATE-master/Matlab/</i>']), ...
                    12, 'center', 'Installing PHATE')
            end
        end
        
        function p=DefineArgs
            p = inputParser;
            defaultMetric = 'euclidean';

            addOptional(p,'csv_file_or_data',[],@(x) ischar(x) || isnumeric(x) || iscell(x));
            addParameter(p,'ndim', 2, @(x) isnumeric(x) && x>=1);
            addParameter(p, 'k', 15, @(x) isnumeric(x) && x>=2);
            addParameter(p, 'a', 40, @(x) isnumeric(x) && x>=0);
            addParameter(p,'t', [], @(x) isempty(x) || ...
                (isnumeric(x) && x>=0)); %interpret 0 as []
            addParameter(p,'t_max', 100, @(x) isnumeric(x) && x>=1);
            addParameter(p, 'npca', 100, @(x) isnumeric(x) && x>=2);
            addParameter(p,'mds_method','mmds', @(x)ismember(x, PhateUtil.MDS_METHODS));
            addParameter(p,'distfun', defaultMetric, @(x)ismember(PhateUtil.METRIC_DICT(x), PhateUtil.METRICS));
            addParameter(p,'distfun_mds', defaultMetric, @(x)ismember(PhateUtil.METRIC_DICT(x), PhateUtil.METRICS));
            addParameter(p,'pot_method', 'log', @(x)ismember(x, PhateUtil.POT_METHODS));
            addParameter(p,'gamma',0.5, @(x) isnumeric(x) && x>=-1 && x<=1);
            addParameter(p,'pot_eps', 1e-7, @(x) isnumeric(x) && x>0 && x<=1);
            addParameter(p,'n_landmarks', 2000, @(x) isempty(x) || (isnumeric(x) && x>= 0));
            addParameter(p,'nsvd', 100, @(x) isnumeric(x) && x>=0);
            addParameter(p,'kernel', []);
            addParameter(p, 'progress_callback', [], ...
                @(x)isequal('function_handle', class(x)));
        end

        function [Y, P, K] = Run(csvFileOrData, varargin)
            PhateUtil.initPaths;

            phateArgs=Args(PhateUtil.DefineArgs);
            args=Args.Str2NumOrLogical(phateArgs.p.Results, varargin);
            for i = 1:length(args)
                if mod(i, 2) == 1 && isKey(PhateUtil.ARG_DICT, args{i})
                    args{i} = PhateUtil.ARG_DICT(args{i});
                    if strcmpi(args{i}, 'distfun')
                        args{i+1} = PhateUtil.METRIC_DICT(args{i+1});
                    end
                end
            end
            args=PhateUtil.Handle0(args);
            [Y, P, K]=phate(csvFileOrData, args{:});
        end

        function initPaths
            pth=fileparts(mfilename('fullpath'));
            pPth=fileparts(pth);
            phatePath=fullfile(pPth, 'phate');
            addpath(phatePath);
        end

        function [argsObj, args, argued, unmatched]...
                =GetArgsWithMetaInfo(varargin)
            if mod(length(varargin),2)==0
                varargin=['sample10k.csv' varargin];
            end
            [args, argued, unmatched, argsObj]=Args.NewKeepUnmatched(...
                PhateUtil.DefineArgs, varargin{:});            
            argsObj.commandPreamble='suh_pipelines';
            argsObj.commandVarArgIn='''pipeline'', ''phate'', ';
            m=mfilename('fullpath');
            argsObj.setSources(@PhateUtil.Run, [m '.m'], m);
            argsObj.setPositionalArgs('csv_file_or_data');
            argsObj.load;
        end
        
        function argsObj=SetArgsMetaInfo(argsObj)
            argsObj.setMetaInfo('k', 'low', 2, 'high', 200, ...
                'is_integer', true, 'label', 'Nearest neighbors', ...
                'text_columns', 2);
            argsObj.setMetaInfo('pot_method', ...
                'type', 'char', ...
                'label', 'PHATE potential distance',...
                'valid_values', PhateUtil.POT_METHODS);
            argsObj.setMetaInfo('mds_method', ...
                'label', 'Method of MDS (multidimensional scaling)',...
                'type', 'char', 'valid_values', ...
                PhateUtil.MDS_METHODS);
            argsObj.setMetaInfo('distfun', ...
                'label', 'Distance function',...
                'type', 'char', 'valid_values', ...
                PhateUtil.METRICS);
            argsObj.setMetaInfo('distfun_mds', ...
                'label', ' MDS distance function',...
                'type', 'char', 'valid_values', ...
                PhateUtil.METRICS);        
            argsObj.setMetaInfo('gamma', 'low', -1, 'high', 1, ...
                'type', 'double', 'text_columns', 5, ...
                'label', 'Gamma potential');
            argsObj.setMetaInfo('pot_eps', 'low', 0, 'high', 1, ...
                'type', 'double', 'text_columns', 8, ...
                'label', 'Epsilon for log potential method');
            argsObj.setMetaInfo('nsvd', 'low', 2, 'high', 300, ...
                'is_integer', true, 'label', '# singular vectors', ...
                'text_columns', 3);
            argsObj.setMetaInfo('a', 'low', 0, 'high', double(intmax), ...
                'label', 'Alpha decaying kernel', 'type', 'double', ...
                'text_columns', 6);
            argsObj.setMetaInfo('t', 'low', 0, 'high', 300, 'label', ...
                '# diffusion steps', 'is_integer', true,...
                'text_columns', 6);
            argsObj.setMetaInfo('t_max', 'low', 50, ...
                'high', intmax, 'label', 'Max diffusion steps', ...
                'is_integer', true, 'text_columns', 3);
            argsObj.setMetaInfo('ndim', 'low', 2, ...
                'high', 100, 'label', '# reduced dimensions', ...
                'is_integer', true, 'text_columns', 2);
            argsObj.setMetaInfo('n_landmarks', 'low', 100, ...
                'high', intmax, 'label', '# landmarks for subsampling', ...
                'is_integer', true, 'text_columns', 2);
            argsObj.setArgGroup({'k', 'n_landmarks', 'distfun', ...
                'ndim', 'pot_method', 'pot_eps', 'mds_method', ...
                'distfun_mds', 't', 't_max', 'a', 'gamma'}, 'Basic settings')
            argsObj.setFileFocus('Unreduced input data', 'csv_file_or_data');
        end
    end
end

