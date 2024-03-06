% MutualInformation: returns mutual information (in bits) of the 'X' and 'Y'
% by Will Dwinnell
%
% I = MutualInformation(X,Y);
%
% I  = calculated mutual information (in bits)
% X  = variable(s) to be analyzed (column vector)
% Y  = variable to be analyzed (column vector)
%
% Note: Multiple variables may be handled jointly as columns in matrix 'X'.
% Note: Requires the 'Entropy' and 'JointEntropy' functions.
%
% Last modified: Nov-12-2006

function I = MutualInformation(X,Y)

if (size(X,2) > 1)  % More than one predictor?
    % Axiom of information theory
    I = JointEntropy(X) + Entropy(Y) - JointEntropy([X Y]);
else
    % Axiom of information theory
    I = Entropy(X) + Entropy(Y) - JointEntropy([X Y]);
end


% God bless Claude Shannon.

% EOF


end


% JointEntropy: Returns joint entropy (in bits) of each column of 'X'
% by Will Dwinnell
%
% H = JointEntropy(X)
%
% H = calculated joint entropy (in bits)
% X = data to be analyzed
%
% Last modified: Aug-29-2006

function H = JointEntropy(X)

    % Sort to get identical records together
    X = sortrows(X);
    
    % Find elemental differences from predecessors
    DeltaRow = (X(2:end,:) ~= X(1:end-1,:));
    
    % Summarize by record
    Delta = [1; any(DeltaRow')'];
    
    % Generate vector symbol indices
    VectorX = cumsum(Delta);
    
    % Calculate entropy the usual way on the vector symbols
    H = Entropy(VectorX);
    
    
    % God bless Claude Shannon.
    
    % EOF
    
    

end

% Entropy: Returns entropy (in bits) of each column of 'X'
% by Will Dwinnell
%
% H = Entropy(X)
%
% H = row vector of calculated entropies (in bits)
% X = data to be analyzed
%
% Example: Measure sample entropy of observations of variables with
%   1, 2, 3 and 4 bits of entropy.
%
% Note: Estimated entropy values are slightly less than true, due to
% finite sample size.
%
% X = ceil(repmat([2 4 8 16],[1e3,1]) .* rand(1e3,4));
% Entropy(X)
%
% Last modified: Nov-12-2006

function H = Entropy(X)

    % Establish size of data
    [n m] = size(X);
    
    % Housekeeping
    H = zeros(1,m);
    
    for Column = 1:m,
        % Assemble observed alphabet
        Alphabet = unique(X(:,Column));
        
        % Housekeeping
        Frequency = zeros(size(Alphabet));
        
        % Calculate sample frequencies
        for symbol = 1:length(Alphabet)
            Frequency(symbol) = sum(X(:,Column) == Alphabet(symbol));
        end
        
        % Calculate sample class probabilities
        P = Frequency / sum(Frequency);
        
        % Calculate entropy in bits
        % Note: floating point underflow is never an issue since we are
        %   dealing only with the observed alphabet
        H(Column) = -sum(P .* log2(P));
    end
    
    
    % God bless Claude Shannon.
    
    % EOF
    
    
end    