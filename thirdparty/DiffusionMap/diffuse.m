function [X, eigenvals, psi, phi] = diffuse(D,eps_val,neigen,t)

% DIFFUSE --- eigendecomposition for diffusion maps
% function [X, eigenvals, psi, phi] = diffuse(D,eps_val,neigen,t)
%
% Input: 
%       D          --- pairwise distances, n-by-n matrix
%       eps_val    --- parameter in weight matrix exp(-D^2/(4*eps_val))
%       neigen     --- number of dimensions of final representation
%                       (default: 95% drop-off in eigenvalue multiplier)
%       t          --- optional time parameter in diffusion map 
%                       (default: model with multiscale geometry)    
%
% Output: 
%       X          --- non-trivial diffusion coords. (entered column-wise)  
%       eigenvals  --- eigenvalues of Markov matrix
%       psi        --- right eigenvectors of Markov matrix
%       phi        --- left eigenvectors
%
%  Ann B Lee, March 2009. Last modified, JWR: 3/23/2009

if (nargin<3), neigen=[]; end % if neigen not defined by user

n=size(D,1);
K = exp(-D.^2/(4*eps_val)); %or equivalently, K=exp(-(D/sigmaK).^2); 
v=sqrt(sum(K)); v=v(:);
A=K./(v*v');   % symmetric graph Laplacian
threshold=5E-6; 
A=sparse(A.*double(A>threshold));  % make matrix sparse to speed up calcs  
if (isempty(neigen))
  [U,S,V]=svds(A,51);  % eigendecomposition of symmetric matrix
  psi=U./(U(:,1)*ones(1,51)); % right eigenv of Markov matrix
  phi=U.*(U(:,1)*ones(1,51)); % left eigenv of Markov matrix
else
  [U,S,V]=svds(A,neigen+1);  % eigendecomposition of symmetric matrix
  psi=U./(U(:,1)*ones(1,neigen+1)); % right eigenv of Markov matrix
  phi=U.*(U(:,1)*ones(1,neigen+1)); % left eigenv of Markov matrix
end
eigenvals=diag(S);

% DIFFUSION COORDINATES
if nargin>3 % fixed scale parameter
  lambda_t=eigenvals(2:end).^t; lambda_t=ones(n,1)*lambda_t';
  if (isempty(neigen)) % use neigen corresponding to 95% drop-off
                         % in lambda_t
    lam = (lambda_t(1,:))/(lambda_t(1,1));
    neigen = find(lam<.05, 1 ); % default number of eigenvalues
    neigen = min([neigen 50]); % use neigen=50 if 95% dropoff not attained
    disp(['Used default value: ',num2str(neigen),' dimensions'])
  end
  X = psi(:,2:neigen+1).*lambda_t(:,1:neigen);  % diffusion coords X
      %= right eigenvectors rescaled with eigenvalues
else  % multiscale geomtry
  lambda_multi = eigenvals(2:end)./(1-eigenvals(2:end));
  lambda_multi=ones(n,1)*lambda_multi';
  if (isempty(neigen)) % use neigen corresponding to
                         % 95% drop-off in lambda_multi
    lam = (lambda_multi(1,:))/(lambda_multi(1,1));
    neigen = find(lam<.05, 1 ); % default number of eigenvalues
    neigen = min([neigen 50]); % use neigen=50 if 95% dropoff not attained
    disp(['Used default value: ',num2str(neigen),' dimensions'])
  end
  X = psi(:,2:neigen+1).*lambda_multi(:,1:neigen);  % diffusion coords X
end

