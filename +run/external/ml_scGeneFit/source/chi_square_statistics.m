function d = chi_square_statistics(XI, XJ)
% Implementation of the Chi^2 distance to use with pdist
% (cf. "The Earth Movers' Distance as a Metric for Image Retrieval",
%      Y. Rubner, C. Tomasi, L.J. Guibas, 2000)
%
% @author: B. Schauerte
% @date:   2009
% @url:    http://cvhci.anthropomatik.kit.edu/~bschauer/

% Copyright 2009 B. Schauerte. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%    1. Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%
%    2. Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the
%       distribution.
%
% THIS SOFTWARE IS PROVIDED BY B. SCHAUERTE ''AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL B. SCHAUERTE OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
% BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
% WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
% OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% The views and conclusions contained in the software and documentation
% are those of the authors and should not be interpreted as representing
% official policies, either expressed or implied, of B. Schauerte.

m = size(XJ, 1); % number of samples of p
p = size(XI, 2); % dimension of samples

assert(p == size(XJ, 2)); % equal dimensions
assert(size(XI, 1) == 1); % pdist requires XI to be a single sample

d = zeros(m, 1); % initialize output array

for i = 1:m
    for j = 1:p
        m = (XI(1, j) + XJ(i, j)) / 2;
        if m ~= 0 % if m == 0, then xi and xj are both 0 ... this way we avoid the problem with (xj - m)^2 / m = (0 - 0)^2 / 0 = 0 / 0 = ?
            d(i, 1) = d(i, 1) + ((XI(1, j) - m)^2 / m); % XJ is the model! makes it possible to determine each "likelihood" that XI was drawn from each of the models in XJ
        end
    end
end