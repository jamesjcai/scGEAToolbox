function [ Kw, ny, Mg, zo,wo] = find_mg( z,lambda,w )
%FIND_MG find the MG other stuff for the prox_l12
%
%   In the litterature, they are mistakes in the formulas. As consquence,
%   this function is a mixed of the formulas found in those 3 papers. In the
%   first paper, the formulas seems correct. The easier to read is most
%   likely the third one.
%
%   Cheers!
%
%   References:
%     M. Kowalski, K. Siedenburg, and M. Dorfler. Social sparsity!
%     neighborhood systems enrich structured shrinkage operators. Signal
%     Processing, IEEE Transactions on, 61(10):2498--2511, 2013.
%     
%     M. Kowalski. Sparse regression using mixed norms. Applied and
%     Computational Harmonic Analysis, 27(3):303--324, 2009.
%     
%     M. Kowalski and B. Torresani. Sparsity and persistence: mixed norms
%     provide simple signal models with dependent coefficients. Signal, image
%     and video processing, 3(3):251--264, 2009.
%     
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/prox/misc/find_mg.php

% Copyright (C) 2012-2016 Nathanael Perraudin.
% This file is part of UNLOCBOX version 1.7.4
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


if nargin<3
    w = ones(size(z));
end


r = abs(z./w);

[Nelg, Ng] = size(r);

[r, ind] = sort(r,1,'descend');

zo = zeros(Nelg,Ng);
wo = zeros(Nelg,Ng);

for jj = 1 : Ng
    zo(:,jj) = z(ind(:,jj),jj);
    wo(:,jj) = w(ind(:,jj),jj);
end

Mg = zeros(Ng,1);
Kw = zeros(Ng,1);
ny = zeros(Ng,1);

for jj = 1 : Ng

    for ii=1:Nelg-1
        temp = lambda * ...
            sum( wo(1:(ii+1),jj).^2 .* ( r(1:(ii+1),jj) - r((ii+1),jj) ) )...
            - r((ii+1),jj);
        if temp >= 0
            Mg(jj) = ii;
            Kw(jj) = sum( wo(1:ii,jj).^2 );
            ny(jj) = norm( zo(1:ii,jj),1);
            break;
        end
    end
    
    % handle exception
    if Mg(jj)==0
        Mg(jj) = 1;
        Kw(jj) = wo(1,jj).^2;
        ny(jj) = abs(zo(1,jj));
    end

end




end

