function [f] = i_waitbar(f)
if nargin < 1 || isempty(f)
    f = waitbar(0, 'Processing your data...');
    pause(.5)
    fprintf('Processing your data...');
    waitbar(.67, f, 'Processing your data...');
    fprintf('... ');
    tic;
    return;
elseif isvalid(f) && strcmp(f.Tag, 'TMWWaitbar')
    toc;
    waitbar(1, f, 'Finishing');
    pause(1);
    % fprintf('.......................done.\n');
    if isvalid(f)
        close(f);
    end
end
