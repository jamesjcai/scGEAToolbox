function cg = unitaryGate(targetQubits, U, NameValueArgs)
%UNITARYGATE Unitary matrix gate
%
%   cg = UNITARYGATE(targetQubits, U) returns a CompositeGate that applies
%   the unitary matrix U to the target qubits up to a global phase.
%
%   cg = UNITARYGATE(___, RotationThreshold=THRES) removes rotation gates
%   with angle magnitude below the threshold THRES. The threshold must be a
%   positive real number or one of the following:
%
%   "auto" - (default) The threshold is 2*pi*eps to remove angles that are
%            approximately 0.
%   "none" - No gates are removed.
%
%   The returned CompositeGate uses about O(4^N) two-qubit gates and many
%   one-qubit gates, where N is the number of target qubits. It is not
%   guaranteed to use the minimal number of gates possible for a given
%   input matrix.
%
%   See also compositeGate, initGate, ucrxGate, ucryGate, ucrzGate

%  Copyright 2023-2024 The MathWorks, Inc.

% References:
% [1] V. V. Shende, S. S. Bullock and I. L. Markov. "Synthesis of quantum-
% logic circuits". IEEE Transactions on Computer-Aided Design of Integrated
% Circuits and Systems. June 2006
% [2] F. Vatan and C. Williams. "Optimal quantum circuits for general
% two-qubit gates". Phys Rev Lett. August 2004
% [3] B. Drury and P. Love. "Constructive quantum Shannon decomposition from
% Cartan involutions". Journal of Physics A. September 2008
arguments
    targetQubits {mustBeVector(targetQubits, 'allow-all-empties'), mustBeInteger, mustBePositive}
    U {mustBeNumeric}
    NameValueArgs.RotationThreshold {quantum.internal.gate.mustBeRotationThreshold(NameValueArgs.RotationThreshold)} = "auto"
end

numQubits = length(targetQubits);
if length(unique(targetQubits))~=numQubits
    error(message("quantum:gates:matchingQubits"))
end

dim = size(U,1);
if ~ismatrix(U) || dim~=size(U,2) || dim~=2^numQubits
    % Input does not represent a valid unitary gate matrix
    error(message("quantum:CompositeGate:invalidUnitarySize", 2^numQubits))
end

% Check input up to a unitary tolerance that scales with the size and datatype
if norm(U'*U - eye(dim), 1) > dim*sqrt(eps('like', U))
    error(message("quantum:CompositeGate:invalidUnitary"))
end

if numQubits==0
    cg = compositeGate(quantum.gate.SimpleGate.empty(0,1), targetQubits, Name="unitary");
    return
end

% Set numeric value for rotation threshold
thres = NameValueArgs.RotationThreshold;
if matlab.internal.math.partialMatch(thres, "none")
    thres = 0;
elseif matlab.internal.math.partialMatch(thres, "auto")
    thres = 2*pi*eps;
end

% Recursively decompose the input matrix
cg = qsdGate(targetQubits, U, thres);
end

%% Gate Decomposition Helpers

function cg = qsdGate(targetQubits, U, thres)
% Recursive Quantum Shannon Decomposition (Theorem 13 [1])

% The type of decomposition depends on the number of qubits
numQubits = length(targetQubits);

if numQubits==1 %#ok<ISCL>
    % Use three 1-qubit rotation gates about the Z and Y axis
    cg = eulerGate(targetQubits, U, thres);
    return
end

if numQubits==2
    % Use several 2-qubit entangling gates between 4 generic 1-qubit rotations
    [cg, isValid] = cartanGate(targetQubits, U, thres);
    if isValid
        return
    end
    % Matrix conditioning for a valid Cartan decomposition failed. Use the
    % generic decomposition below which uses an additional 2 two-qubit
    % gates.
end

% Decompose into 4 smaller unitary matrices on 1 less qubit, and alternate
% them between 3 uniform rotation gates. These act on the top qubit and are
% controlled on the bottom qubits.
trgt = 1;
ctrls = 2:numQubits;

% Use the Cosine-Sine decomposition to factor U into the form:
%
% U = blkdiag(V1, V2)*[C -S; S C]*blkdiag(W1, W2)
%
% The V and W are smaller unitary matrices, and C and S are diagonal matrices
% used to compute the rotation angles for the uniform rotation about the Y axis.
[V1, V2, W1, W2, ucryAngles] = csd2by2(U);

% Construct the uniform Y-axis rotation with czGates instead of cxGates
ucryWithCZ = quantum.internal.gate.ucrgates(ucryAngles, @ryGate, @czGate, thres);
if ~isempty(ucryWithCZ) && ...
        (ucryWithCZ(end).Type=="cz" && ucryWithCZ(end).ControlQubits==1)
    % The czGate is diagonal, so it can be merged with the generic multiplexed
    % unitary gate on the right. Therefore, 1 czGate can be removed at each
    % recursive step by applying it's tensored matrix to blkdiag(V1, V2).
    % Only the czGate initially controlled on the first qubit is merged.
    % This only requires multiplication to the columns of V2 since the
    % top-left block applied to V1 is the identity.
    ucryWithCZ(end) = [];
    V2 = V2 .* repelem([1 -1], 2^(numQubits-2));
end

ucryWithCZ = compositeGate(ucryWithCZ, [ctrls trgt], Name="ucry");

% Factor these smaller unitary matrices into unitary matrices G using
%
% blkdiag(U1,U2) = blkdiag(G2,G2)*blkdiag(D,D')*blkdiag(G1,G1)
%
% D is diagonal and used to compute rotation angles for the uniform
% rotation about the Z axis.
[G2, G1, ucrzAnglesL] = demultiplex(W1, W2);
[G4, G3, ucrzAnglesR] = demultiplex(V1, V2);

% Decompose the smaller unitary matrices to construct the
% final circuit that applies the input unitary matrix.
gates = [
    qsdGate(ctrls, G1, thres)
    ucrzGate(ctrls, trgt, ucrzAnglesL, RotationThreshold=thres)
    qsdGate(ctrls, G2, thres)
    ucryWithCZ
    qsdGate(ctrls, G3, thres)
    ucrzGate(ctrls, trgt, ucrzAnglesR, RotationThreshold=thres)
    qsdGate(ctrls, G4, thres)
    ];
cg = compositeGate(gates, targetQubits, Name="unitary");
end

function cg = eulerGate(targetQubit, U, thres)
% Euler rotations in the basis ZYZ (Theorem 1 [1])
% U = exp(i*phi)*RZ(alpha)*RY(beta)*RZ(gamma)

% Scale U to have determinant 1
U = U / sqrt(det(U));

gamma = angle(U(2,2))-angle(U(2,1));
beta = 2*atan2(abs(U(2,1)), abs(U(1,1)));
alpha = angle(U(2,2))+angle(U(2,1));

gates = [rzGate(1, gamma) ryGate(1, beta) rzGate(1, alpha)];
gates = quantum.internal.gate.filterRotationGates(gates, thres);
cg = compositeGate(gates, targetQubit, Name="unitary");
end


function [cg, success] = cartanGate(targetQubits, U, thres)
% KAK Decomposition (Section 3.2 [3])
% U = K1*A*K2 = kron(Va,Vb)*B*D*B'*kron(Ua,Ub)

% Scale U to have determinant 1
U = U / detUsingQR(U)^(1/4);

% Transform U to the magic basis
B = 1/sqrt(2).*[1 1j 0 0; 0 0 1j 1; 0 0 1j -1; 1 -1j 0 0];
Up = B'*U*B;

M2 = Up.'*Up;
[P, D2] = schur(M2, 'complex');
d2 = diag(D2);

% Scale P to be real by rotating each column by conjugate of its max angle.
% In practice, this scales P to be a real orthogonal matrix but is not
% guaranteed, so check if this was successful in the next step.
for jj = 1:size(P,2)
    [~, idx] = max(abs(P(:,jj)));
    P(:,jj) = conj(sign(P(idx,jj))) * P(:,jj);
end
P(:,1) = P(:,1) / detUsingQR(P);

success = norm(imag(P), 1) < 1e-14;
if ~success
    % Decomposition was not successful because imaginary components exceed
    % tolerance.
    cg = [];
    return
end
P = real(P);

% The square root results in D having determinant -1 or 1. But it must be
% 1 since Kp computed below must also have determinant 1 for K1 to be
% expressed as a Kronecker product.
d = sqrt(d2);
if real(prod(d)) < 0
    % Make the determinant 1 by swapping the sign of an element since the
    % determinant of diagonal matrix is the product of it elements.
    d(1) = -d(1);
end

% Solve system of angles and construct the entangling circuit from
% Figure 6 [2]
phases = angle(d);
coeffs = [1 -1 1 1; -1 1 1 1; 1 1 -1 1; -1 -1 -1 1];
params = coeffs\phases;

entangleGate = [
    rzGate(2, -pi/2)
    cxGate(2, 1)
    rzGate(1, mod(-2*params(3)+(pi/2), 2*pi))
    ryGate(2, mod(-(pi/2)+2*params(1), 2*pi))
    cxGate(1, 2)
    ryGate(2, mod(-2*params(2)+(pi/2), 2*pi))
    cxGate(2, 1)
    rzGate(1, pi/2)
    ];
entangleGate = quantum.internal.gate.filterRotationGates(entangleGate, thres);

% Undo magic basis transformation so K1 and K2 can be expressed as a
% Kronecker product of two unitary matrices
Dinv = diag(1./d);
Kp = Up*P*Dinv*P';
K1 = B*Kp*P*B';
K2 = B*P'*B';

[Va, Vb] = nearestKron(K1);
[Ua, Ub] = nearestKron(K2);

% Final circuit is the entangling gates between the 4 generic 1-qubit
% rotations shown in Figure 1 [3].
gates = [
    eulerGate(1,Ua, thres)
    eulerGate(2,Ub, thres)
    entangleGate
    eulerGate(1,Va, thres)
    eulerGate(2,Vb, thres)
    ];
cg = compositeGate(gates, targetQubits, Name="unitary");
end

%% Matrix Decomposition Helpers

function [B,C] = nearestKron(A)
% Nearest Kronecker Decomposition
% A = kron(B,C)

% Construct permuted version of A
pA = reshape(A, [2 2 2 2]);
pA = permute(pA, [2 4 1 3]);
pA = reshape(pA, [4 4]);

% Rank-1 approximation
[U,sv,V] = svd(pA);
vecB = sqrt(sv(1,1))*U(:,1);
vecC = sqrt(sv(1,1))*conj(V(:,1));

B = reshape(vecB, [2 2]);
C = reshape(vecC, [2 2]);
end


function [V1, V2, W1, W2, angles] = csd2by2(U)
% 2-by-2 Cosine-Sine Decomposition (Theorem 10 [1])
% U = blkdiag(V1, V2)*[C -S; S C]*blkdiag(W1, W2)

sz = size(U,1);
half = sz/2;

[V1, V2, W1, C, S] = gsvd(U(1:half, 1:half), U(half+1:end, 1:half));
W1 = W1';
W2 =  [-V1*S; V2*C]'*U(:, half+1:end);

angles = 2*atan2(diag(S), diag(C));
end


function [V, W, angles] = demultiplex(U1, U2)
% Multiplexor Decomposition (Theorem 12 [1])
% blkdiag(U1, U2) = blkdiag(V, V)*blkdiag(D, D')*blkdiag(W, W)

[V, D2] = schur(U1*U2', 'complex');

d = sqrt(diag(D2));
W = d.*(V'*U2);

angles = -2*angle(d);
end


function d = detUsingQR(A)
% Compute the determinant of input matrix A by way of the QR decomposition.
% This gives a numerical safe result when A is an orthogonal matrix, where
% the LU decomposition used in DET may return wrong result.
%
% For an example, det(orth(randn(5e3,'like',1i))) returns 0 in most cases,
% when the correct result (which detUsingQR returns) has absolute value 1.
%
% This function also works for non-orthogonal A, but doesn't have any clear
% advantages over det in those cases.

% Compute QR decomposition of A, storing the LAPACK output format with
% combined R and Householder reflectors, and where tau is the scalar
% multiplier for each Householder reflector.
% A = Q * R, where R = triu(H) and Q = Q1 * ... * Qn,
%                  with Qi = I - tau(i) * vi * vi'
% v_i can be retrieved from H, but we don't need it here, since we only
% care about det(Qi).
[H, tau] = matlab.internal.decomposition.compactQR(A);

% det(Qi) = det(I - tau(i) * vi * vi') = 1 - tau(i)*vi'*vi
%
%  We know abs(det(Qi)) == 1, and can deduce vi'*vi from this:
%      vi'*vi = 2*real(tau(i)) / abs(tau(i))^2
%
% Inserting this above, we find
% det(Qi) = 1 - 2*tau(i)*real(tau(i)) / abs(tau(i))^2
%         = 1 - 2 * real(s(i)) * s(i), with s(i) = sign(tau(i)).
%
% If abs(tau(i)) ~= 0, this could be simplified further to -s(i)^2, but in
% cases where no reflection is needed, tau(i) is set to zero.
s = sign(tau);
detOfHouseholderReflectors = 1 - 2*real(s).*s;
diagR = diag(H);

% Compute det(A) = det(Q) * det(R) = det(Q1) *...* det(Qn) * prod(diag(R))
% Both vectors above have all elements with absolute value 1, but with
% varying angles. Computing their product gives us the expected determinant
% in a safe way.
d = prod(detOfHouseholderReflectors) * prod(diagR);
end