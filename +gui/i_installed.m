function y = i_installed(name, parentfig)

if nargin<2, parentfig = []; end

y = false;
switch name
    case 'stats'
    if (ispc && (~license('test','statistics_toolbox') || isempty(which('grp2idx.m')))) || ...
       (~ispc && ~any(strcmp('Statistics and Machine Learning Toolbox', {ver().Name})))
        
        gui.myWarndlg(parentfig, ['SCGEATOOL requires Statistics ' ...
            'and Machine Learning Toolbox.'], ...
            'Missing Dependencies');
        if strcmp(gui.myQuestdlg(parentfig, 'Learn how to install Statistics and Machine Learning Toolbox?'),'Yes')
            web('https://www.mathworks.com/help/matlab/matlab_env/get-add-ons.html');
            web('https://www.mathworks.com/videos/add-on-explorer-106745.html');
        end
    else
        y = true;
    end
    otherwise
        error('Unknow tag.');
end