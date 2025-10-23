function [cs, tflist, gcommon, numtargetgenes] = sc_tfactivity(X, g, ...
    Ttfgn, speciestag, methodid)
% The activity level of a transcription factor (TF) in a given cell is the
% extent to which it is exerting its regulatory potential on its target
% genes.
% https://academic.oup.com/bioinformatics/article/37/9/1234/5949002
%
% [cs,tflist]=sc_tfactivity(X,g);
% CS - is an m-by-n matrix of activity scores for m TFs and n cells.
% TFlist - list of TF genes

if nargin < 2, error('USAGE: [cs, tflist] = sc_tfactivity(X, g);'); end
if nargin < 5 || isempty(methodid), methodid = 4; end
if nargin < 4 || isempty(speciestag), speciestag = 'hs'; end
% if nargin < 3 || isempty(Ttfgn) % tf-by-gene matrix T from database
%     %folder=fileparts(mfilename('fullpath'));
%     %wrkpth=fullfile(folder,'assets',filesep,'DoRothEA_TF_Target_DB',filesep);
%     pw1 = fileparts(mfilename('fullpath'));
%     switch lower(speciestag)
%         case {'hs', 'human'}
%             %fname=[wrkpth 'dorothea_hs.mat'];
%             fname = fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat');
%         case {'mm', 'mouse'}
%             %fname=[wrkpth 'dorothea_mm.mat'];
%             fname = fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_mm.mat');
%         otherwise
%             error('TF database is not available for the species.');
% 
%     end
%         fprintf('\nReading ... %s.\n', fname);
%             load(fname, 'T');
%             Ttfgn = T(T.mor > 0, :);
%             fprintf('Only positive regulatory relationships are used.\n');
% end

if nargin < 3 || isempty(Ttfgn)
    pw1 = fileparts(mfilename('fullpath'));
    if strcmpi(speciestag, 'hs') || strcmpi(speciestag, 'human')
        fname = fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat');
    elseif strcmpi(speciestag, 'mm') || strcmpi(speciestag, 'mouse')
        fname = fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_mm.mat');
    else
        error('TF database is not available for the species.');
    end
    
     if exist(fname, 'file')
         fprintf('\nReading ... %s.\n', fname);
         load(fname, 'T');
         if ~exist('T', 'var')
             error('Variable T not found in the file %s.', fname);
         end
     else
         error('File %s does not exist.', fname);
     end
    Ttfgn = T(T.mor > 0, :); % Filter positive regulatory relationships
    fprintf('Only positive regulatory relationships are used.\n');
end


if ~isnumeric(X) || ~ismatrix(X)
    error('X must be a numeric matrix.');
end
if ~iscellstr(g) && ~isstring(g)
    error('g must be a cell array of strings.');
end
% if ~all(isfield(Ttfgn.Properties.VariableNames, {'tf', 'target', 'mor'}))
if ~all(contains({'tf', 'target', 'mor'}, Ttfgn.Properties.VariableNames))
    error('Ttfgn must contain the fields: tf, target, and mor.');
end


    try
        if issparse(X), X = full(X); end
    catch
        warning('Keep using sparse X.');
    end

    if methodid ~= 1 % method 1 UCell is rank-based, normalization is unnecessary
        [X] = sc_norm(X);
        [X] = log1p(X);
    end

    [gid, gnlist] = findgroups(string(Ttfgn.target));
    [tid, tflist] = findgroups(string(Ttfgn.tf));
    t = zeros(max(tid), max(gid));
    t(sub2ind([max(tid), max(gid)], tid, gid)) = Ttfgn.mor;
    

    fprintf('Using the Dorothea dadtabase that contains %d TFs and %d targets.\n', ...
        size(t, 1), size(t, 2));

    % size(t)
    % assignin('base','t2',t);

    %t2=zeros(max(tid),max(gid));
    %for k=1:length(gid)
    %    t2(tid(k),gid(k))=T.mor(k);
    %end
    %isequal(t,t2)
    % T=T(T.mor>0,:);    % only consider positive regulation
    %[t]=crosstab(T.tf,T.target);   % TF-by-target regulagory relationship
    %matrix if only positive regulation

    [~, k, l] = intersect(upper(g), upper(gnlist));
    t = t(:, l); % tf-by-gene matrix
    X = X(k, :); % gene-by-cell matrix
    fprintf('Using %d target genes that are present in the data.\n', size(t, 2));

    if nargout > 2, gcommon = g(k); end

    switch methodid
        case 1 % UCell method  see also: sc_cellscore_ucell

            cs = zeros(size(t, 1), size(X, 2));
            R = tiedrank(-X);
            R(R > 1500) = 1500 + 1;
            numtargetgenes = zeros(size(t, 1), 1);
            for k = 1:size(t, 1)
                idx1 = t(k, :) > 0;
                n1 = sum(idx1);
                if n1 > 0
                    u = sum(R(idx1, :)) - (n1 * (n1 - 1)) / 2;
                    cs(k, :) = 1 - u / (n1 * 1500);
                    numtargetgenes(k) = n1;
                end
            end
            cs(cs < 0) = 0;

            % idx1 = t > 0; % Logical matrix for positive relationships
            % numtargetgenes = sum(idx1, 2); % Precompute target gene counts
            % R = tiedrank(-X);
            % R(R > maxRank) = maxRank + 1;
            % u = sum(idx1 .* R, 2) - (numtargetgenes .* (numtargetgenes - 1)) / 2;
            % cs = 1 - u ./ (numtargetgenes * maxRank);
            % cs(cs < 0) = 0;

        case 2 % matrix multiplication method
            cs = t * X;
            numtargetgenes = sum(t > 0, 2);
        case 3 % matrix inverse method
            cs = pinv(t') * X;
            numtargetgenes = sum(t > 0, 2);
        case 4 % nnmf method
            disp('ref: Bioinformatics, Volume 37, Issue 9, 1 May 2021, Pages 1234–1245,');
            disp('https://doi.org/10.1093/bioinformatics/btaa947');
            disp('PMID: 33135076 PMCID: PMC8189679 DOI: 10.1093/bioinformatics/btaa947');
            n = size(t, 1);
            v.WRfixed = n;
            v.W = t.';
            [~, cs] = NMF(X, n, 100, 0, v);
            numtargetgenes = sum(t > 0, 2);
    end
end


function [W, H, bDsave] = NMF(V, R, Niter, beta, initialV)
    % [W,H, bDsave] = NMF(V,R,Niter,beta,initialV)
    %    NMF with beta divergence cost function.
    %Input :
    %   - V : power spectrogram to factorize (a MxN matrix)
    %   - R : number of templates
    %   - Niter : number of iterations
    %   - beta (optional): beta used for beta-divergence (default : beta = 0, IS divergence)
    %   - initialV (optional) : initial values of W, H (a struct with
    %   fields W and H)
    %Output :
    %   - W : frequency templates (MxR array)
    %   - H : temporal activation
    %   - bDsave : evolution of beta divergence
    %
    % Copyright (C) 2010 Romain Hennequin

    % https://github.com/romi1502/NMF-matlab
    disp('NMF implementation: https://github.com/romi1502/NMF-matlab');

    verbose = false;

    eta = 1;

    % size of input spectrogram
    M = size(V, 1);
    N = size(V, 2);

    % initialization
    if nargin == 5
        if isfield(initialV, 'H')
            H = initialV.H;
        else
            H = rand(R, N);
        end
        if isfield(initialV, 'W')
            W = initialV.W;
        else
            W = rand(M, R);
        end

        if isfield(initialV, 'HRfixed')
            HRfixed = initialV.HRfixed;
        else
            HRfixed = 0;
        end

        if isfield(initialV, 'WRfixed')
            WRfixed = initialV.WRfixed;
        else
            WRfixed = 0;
        end


    else
        H = rand(R, N);
        W = rand(M, R);
        HRfixed = 0;
        WRfixed = 0;

        if nargin == 3
            beta = 0;
        end
    end

    % array to save the value of the beta-divergence
    bDsave = zeros(Niter, 1);

    % computation of Lambda (estimate of V) and of filters repsonse
    Lambda = W * H;

    % Waitbar
    message = ['Computing NMF .... iteration : 0/', int2str(Niter), ' completed'];
    h = waitbar(0, message);


    % iterative computation
    for iter = 1:Niter

        %     % compute beta divergence and plot its evolution
        bDsave(iter) = betaDiv(V+eps, Lambda+eps, beta);

        % update of W
        if not(WRfixed)
            W = W .* ((Lambda.^(beta - 2) .* V) * H' + eps) ./ ((Lambda.^(beta - 1)) * H' + eps);
        else
            W(:, WRfixed+1:end) = W(:, WRfixed+1:end) .* ((Lambda.^(beta - 2) .* V) * H(WRfixed+1:end, :)' + eps) ./ ((Lambda.^(beta - 1)) * H(WRfixed+1:end, :)' + eps);
        end

        % recomputation of Lambda (estimate of V)
        Lambda = W * H + eps;


        % update of H
        if not(HRfixed)
            H = H .* (W' * (Lambda.^(beta - 2) .* V) + eps) ./ (W' * (Lambda.^(beta - 1)) + eps);
        else
            H(1:HRfixed, :) = H(1:HRfixed, :) .* (W(:, 1:HRfixed)' * (Lambda.^(beta - 2) .* V) + eps) ./ (W(:, 1:HRfixed)' * (Lambda.^(beta - 1)) + eps);
        end
        % recomputation of Lambda (estimate of V)
        Lambda = W * H + eps;


        message = ['Computing NMF. iteration : ', int2str(iter), '/', int2str(Niter)];
        if verbose
            disp(message);
        end
        waitbar(iter/Niter, h, message);
    end

    % % normalization
    % for r0=1:R
    %     % normalization of templates
    %     chosenNorm = 2;
    %     normW = norm(W(:,r0),chosenNorm);
    %     H(r0,:) = normW*H(r0,:);
    %     W(:,r0) = W(:,r0)/normW;
    % end

    close(h)
    % close
end


function bD = betaDiv(V, Vh, beta)
    if beta == 0
        bD = sum((V(:) ./ Vh(:))-log(V(:)./Vh(:))-1);
    elseif beta == 1
        bD = sum(V(:).*(log(V(:)) - log(Vh(:)))+Vh(:)-V(:));
    else
        bD = sum(max(1/(beta * (beta - 1))*(V(:).^beta + (beta - 1) * Vh(:).^beta - beta * V(:) .* Vh(:).^(beta - 1)), 0));
    end
end
