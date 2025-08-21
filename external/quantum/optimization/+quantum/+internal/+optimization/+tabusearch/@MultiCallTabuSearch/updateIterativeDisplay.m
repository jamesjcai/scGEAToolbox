function updateIterativeDisplay(obj, numTabuCalls, bestFval, elapsedTime, tabuIter)
%UPDATEITERATIVEDISPLAY Update iterative display
%
%   UPDATEITERATIVEDISPLAY(OBJ) updates the iterative display if requested.

if ~strcmp(obj.Display, 'iter')
    return
end

% Display header if needed. 
if rem(numTabuCalls, 30) == 0 
    fprintf('\nTabu call    Best fval   Stall time   Tabu iterations\n');    
end

% Display entry
fprintf("%9d   %10.4g   %10.4g   %9d\n", numTabuCalls, bestFval, elapsedTime, tabuIter);

end