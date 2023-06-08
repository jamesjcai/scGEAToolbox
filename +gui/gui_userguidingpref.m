function [Button]=gui_userguidingpref()


mfolder = fileparts(mfilename('fullpath'));
fullfile(mfolder,'..', 'resources', 'userguiding.png')
[imgx, mapx] = imread(fullfile(mfolder,'..', 'resources', 'userguiding.gif'));
ptImagex = ind2rgb(imgx, mapx);
%ptx.CData = ptImagex;

S.Default = 'Yes';
S.IconString = 'custom';
S.IconData = ptImagex;
Button=buttondlg('User Onboarding Toolbar helps you onboard with walkthroughs to prompt the right in-app experience. Show User Onbarding Toolbar again?', ...
    '','Yes','No','Cancel',S);
