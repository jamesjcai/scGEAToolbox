function [Button]=gui_userguidingpref(onlyinfo)

if nargin<1, onlyinfo=true; end

mfolder = fileparts(mfilename('fullpath'));
%fullfile(mfolder,'..', 'resources', 'userguiding.png')

try
    [imgx, mapx] = imread(fullfile(mfolder,'..', 'resources', 'userguiding.gif'));
    ptImagex = ind2rgb(imgx, mapx);
    S.IconData = ptImagex;
catch
    S.IconData = uint8(rand(50,50,3).*256);
end
if ~ispref('scgeatoolbox','useronboardingtoolbar')
   setpref('scgeatoolbox','useronboardingtoolbar',true);
end
if onlyinfo
    S.Default = 'OK';
    S.IconString = 'custom';
    Button=buttondlg(['User Onboarding Toolbar helps you onboard ' ...
        'with walkthroughs to prompt the right in-app experience.'], ...
        'Thanks for choosing SCGEATOOL','OK',S);
else
    S.Default = 'Yes';
    S.IconString = 'custom';
    Button=buttondlg(['User Onboarding Toolbar helps you onboard ' ...
        'with walkthroughs to prompt the right in-app experience. ' ...
        'Show User Onbarding Toolbar next time?'], ...
        'Thanks for choosing SCGEATOOL','Yes','No','Cancel',S);
    switch Button
        case 'Yes'
            setpref('scgeatoolbox','useronboardingtoolbar',true);
        case 'No'
            setpref('scgeatoolbox','useronboardingtoolbar',false);
        case 'Cancel'
    end
end
