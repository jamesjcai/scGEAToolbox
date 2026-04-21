function [tf, msg] = i_isreportgenavailable(feature)
%I_ISREPORTGENAVAILABLE Check whether MATLAB Report Generator is available.

if nargin < 1 || isempty(feature)
    feature = 'dom';
end

feature = lower(string(feature));
tf = false;
msg = 'MATLAB Report Generator is not installed or licensed.';

hasLicense = license('test', 'MATLAB_Report_Gen') || ...
    license('test', 'matlab_report_gen');

switch feature
    case "ppt"
        hasClasses = ~isempty(which('mlreportgen.ppt.Presentation'));
    otherwise
        hasClasses = ~isempty(which('mlreportgen.dom.Document'));
end

tf = hasLicense && hasClasses;

if tf
    msg = '';
elseif ~hasClasses
    msg = 'MATLAB Report Generator is not installed.';
elseif ~hasLicense
    msg = 'MATLAB Report Generator is installed, but no license is available.';
end

end
