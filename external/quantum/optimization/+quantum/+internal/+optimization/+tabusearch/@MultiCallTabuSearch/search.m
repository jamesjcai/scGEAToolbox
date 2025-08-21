function [bestx, bestf, exitflag, output] = search(obj, C, constantForDisplay)
%SEARCH Search for solutions to a QUBO
%
%   [BESTX, BESTF, EXITFLAG, OUTPUT] = SEARCH(OBJ, C, CONSTANTFORDISPLAY)
%   searches for solutions to the supplied QUBO by calling tabu search
%   multiple times. The iterative display reflects the problem constant via
%   the addition of CONSTANTFORDISPLAY.
%
%   Note: BESTF should not include the CONSTANTFORDISPLAY. This is added
%   back in the ABSTRACTALGORITHM superclass.

%   Copyright 2022-2023 The MathWorks, Inc.


% QuboProblem minimizes, whereas TabuSearch maximizes
C = -C;

% Start time
startTime = tic;

% Get the initial point
x = obj.initialPoint(size(C,1));

% Steps 0-2: Prepare and transform the QUBO
[Cp,dp,x,bestx,bestf] = prepareQuboForTabu(obj, C, x);

% Prepare for cpp driver if used
if obj.UseCppTabuEngine
    fixedStruct = createFixedStruct(Cp);
end

% Time first best
timeFirstBest = 0;

% Record number of times tabu search is called
numTabuSearchCalls = 0;

% Record number of iterations from each tabu search call
allTabuIterations = [];

% Record solutions and function values returned from each tabu search call
allX = bestx;
allFunctionValues = bestf;

% Max stall time
if isempty(obj.MaxStallTime)
    numVars = size(C, 1);
    if numVars <= 250
        maxStallTime = 0.5;
    elseif numVars <= 500
        maxStallTime = 2;
    elseif numVars <= 1000
        maxStallTime = 3;
    elseif numVars <= 2500
        maxStallTime = 3 + (7/1500)*(numVars - 1000);
    else
        maxStallTime = 10;
    end
else
    maxStallTime = obj.MaxStallTime;
end

% Initialize stopState
stopState.LongestStallTime = 0;
stopState.MaxTimeExceeded = false;
stopState.MaxStallTimeExceeded = false;
stopState.NumTabuCallsMet = false;
stopState.Done = false;

% Set up iterative display
updateIterativeDisplay(obj, numTabuSearchCalls, -bestf + constantForDisplay, 0, 0);


numSingleTabuCalls = 0;

elapsedTime = toc(startTime);

stopState = checkStopState(obj, stopState, elapsedTime, timeFirstBest, maxStallTime, numSingleTabuCalls);
while ~stopState.Done

    % Step 3: Call tabusearch
    initf = x'*C*x;
    oldbestf = bestf;
    if obj.UseCppTabuEngine
        [bestx, bestf, Cp, dp, currx, iter] = callTabuOnTransformedProblemCpp(obj, Cp, dp, x, initf, bestx, bestf, fixedStruct, elapsedTime);
    else
        [bestx, bestf, Cp, dp, currx, iter] = callTabuOnTransformedProblem(obj, Cp, dp, x, initf, bestx, bestf, startTime);
    end
    numSingleTabuCalls = numSingleTabuCalls + 1;
    allTabuIterations = [allTabuIterations, iter]; %#ok<AGROW>
    allX = [allX, bestx]; %#ok<AGROW>
    allFunctionValues = [allFunctionValues, bestf]; %#ok<AGROW>
    if obj.ShowDiagnostics
        disp("fcalcBest = " + (bestx'*C*bestx) + " : fvalBest = " + bestf)
    end
    numTabuSearchCalls = numTabuSearchCalls + 1;
    elapsedTime = toc(startTime);
    if bestf > oldbestf
        timeFirstBest = elapsedTime;
    end

    % Step 4.  Check if a stopping criterion is satisfied.
    stopState = checkStopState(obj, stopState, elapsedTime, timeFirstBest, maxStallTime, numSingleTabuCalls);

    % Step 5. Generate state for next iteration. Show iterative display if
    % requested. 
    if ~stopState.Done
        [Cp, dp, x] = generateNextCallState(obj, Cp, dp, currx);
        % Calling generateNextCallState could mean that MaxTime or
        % MaxStallTime are exceeded. Call checkStopState again to make sure
        % time limits aren't exceeded.
        elapsedTime = toc(startTime);
        stopState = checkStopState(obj, stopState, elapsedTime, timeFirstBest, maxStallTime, numSingleTabuCalls);
    end
    
    % Negate function values as quboProblem minimizes
    updateIterativeDisplay(obj, numTabuSearchCalls, -bestf + constantForDisplay, stopState.StallTime, iter);

end

% TabuSearch maximizes, so we need to negate the function values
bestf = -bestf;
allFunctionValues = -allFunctionValues;

% Remove duplicates from allX/allFunctionValues
[~, idx] = unique(allX', 'rows');
allX = allX(:, idx);
allFunctionValues = allFunctionValues(idx);

% Prepare exitflag and output structure
exitflag = 0;
output.ElapsedTime = elapsedTime;
output.TimeFirstBest = timeFirstBest;
output.TabuIterations = allTabuIterations;
output.NumTabuSearchCalls = numTabuSearchCalls;
output.AllX = allX;
output.AllFunctionValues = allFunctionValues;
output.LongestStallTime = stopState.LongestStallTime;
if obj.ReturnTransformedProblem
    % Intended for internal testing and debugging only.
    output.Cp = Cp;
    output.dp = dp;
    output.currx = currx;
end
if stopState.MaxTimeExceeded    
    output.Message = getString(message('quantum:annealing:TabuSearch:ExitMaxTimeExceeded'));
elseif stopState.MaxStallTimeExceeded
    output.Message = getString(message('quantum:annealing:TabuSearch:ExitMaxStallTimeExceeded'));
else
    output.Message = "";
end

% Show final message
if ~strcmp(obj.Display, 'off') 
    fprintf("\n%s\n", output.Message);
end

end

function stopState = checkStopState(obj, stopState, elapsedTime, timeFirstBest, maxStallTime, numSingleTabuCalls)

stopState.StallTime = elapsedTime - timeFirstBest;
if stopState.StallTime > stopState.LongestStallTime
    stopState.LongestStallTime = stopState.StallTime;
end
stopState.MaxTimeExceeded = elapsedTime > obj.MaxTime;
stopState.MaxStallTimeExceeded = stopState.StallTime > maxStallTime;
stopState.NumTabuCallsMet = numSingleTabuCalls == obj.MaxSingleTabuCalls;
stopState.Done = stopState.MaxTimeExceeded || stopState.MaxStallTimeExceeded || stopState.NumTabuCallsMet;

end
