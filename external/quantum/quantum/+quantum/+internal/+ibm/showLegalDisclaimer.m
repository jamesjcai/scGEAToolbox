function showLegalDisclaimer()

%   Copyright 2023 The MathWorks, Inc.

persistent wasDisplayed;
mlock

if isempty(wasDisplayed)
    % Display once per new MATLAB session
    msg = ...
        "The availability of IBM quantum computing systems for submission of quantum inputs "+newline+...
        "by you (the End User) is provided AS IS and is limited to experimental research, "+newline+...
        "education, testing, evaluation and feedback purposes by End Users, and no commercial or "+newline+...
        "production use by End Users is permitted. Results are not guaranteed.";
    disp(msg)
    wasDisplayed = true;
end

