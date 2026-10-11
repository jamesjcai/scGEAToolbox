function [A0, XM0denoised] = i_td1(XM0, methodid)
% TD - tensor decomposition for denoising
%
% input:  XM0 - k multi-layer network array (n x n x k)
%         methodid - 1 for CP (default), 2 for Tucker
% output: A0 - n x n adjacency matrix of denoised network
%         XM0denoised - the denoised n x n x k stack
%
% With one output, A0 comes straight from the factors - the mean of the
% estimated slices - so the n x n x k estimate (720 MB at 3000 genes and 10
% networks) is never built. Neither decomposition needs Tensor Toolbox; see
% ten.i_cpalstt and ten.i_tuckeralstt.
if nargin < 2, methodid = 1; end
if nargout < 2
    switch methodid
        case 1
            [lambda, U] = ten.i_cpalstt(XM0, size(XM0, 3));
            A0 = (U{1}.*(lambda.'.*mean(U{3}, 1)))*U{2}.';
        case 2
            [core, U] = ten.i_tuckeralstt(XM0, 5);
            A0 = U{1}*sum(core.*reshape(mean(U{3}, 1), 1, 1, []), 3)*U{2}.';
        otherwise
            error('ten:i_td1:badMethod', 'METHODID must be 1 (CP) or 2 (Tucker).');
    end
    return
end
switch methodid
    case 1
        Xhat0 = ten.do_td_cp(XM0);
    case 2
        Xhat0 = ten.do_td_tucker(XM0);
    otherwise
        error('ten:i_td1:badMethod', 'METHODID must be 1 (CP) or 2 (Tucker).');
end
A0 = mean(Xhat0, 3);
XM0denoised = Xhat0;
end
