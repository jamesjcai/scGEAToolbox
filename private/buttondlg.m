function ButtonName=buttondlg(Question,Title,Btn1,Btn2,Btn3,Default)
%BUTTONDLG Button dialog box with variable icon.
% *************************************************************************
% ******** Modified function help *****************************************
% *************************************************************************
% buttondlg is a simple modification of MATLABs function QUESTDLG which
% enables the user to change the icon in the dialog box, known from the 
% MATLAB function MSGBOX. Available icons are the same known from MSGBOX
% ('none','warn','error','help','custom'), but also 'quest' (known from the
% original questdlg) and 'bulb'. If you choose 'custom', you also have to
% specify IconData and optional IconCMap, also known form fhe function
% MSGBOX. Best icon size is [50,50]. To submit icon information use the
% last input variable 'Default', similar to 'OPTIONS' in QUESTDLG. This
% variable needs to be a struct with some or all of the following fields:
% Default.Interpreter --> known from QUESTDLG, default is 'none'
% Default.Default --> known from QUESTDLG, default is ''
% Default.IconString --> icon type known from MSGBOX, default is 'quest'
% Default.IconData --> icon image data known from MSGBOX
% Default.IconCMap --> icon color map known from MSGBOX
% Example 1:
% Button=buttondlg('System Error!','Error','Stop','Continue', ...
%       struct('Default','Continue','IconString','error'));
% Example 2:
% Button=buttondlg('Do you like the color icon?','Question','Yes',...
%       'Maybe','No',struct('Default','Yes','IconString','custom',...
%       'IconData',uint8(rand(50,50,3).*256)));
% % is the same as:
% S.Default = 'Yes';
% S.IconString = 'custom';
% S.IconData = uint8(rand(50,50,3).*256);
% Button=buttondlg('Do you like the color icon?','Question','Yes',...
%       'Maybe','No',S)
% Example 3:
% options.Interpreter = 'tex';
% options.Default = 'Don''t know';
% options.IconString = 'none';
% qstring = 'Is \Sigma(\alpha - \beta) < 0?';
% title = 'Boundary Condition';
% Button = buttondlg(qstring,title,'Yes','No','Don''t know',options)
% *************************************************************************
% ******** Original function help *****************************************
% *************************************************************************
%  ButtonName = QUESTDLG(Question) creates a modal dialog box that
%  automatically wraps the cell array or string (vector or matrix)
%  Question to fit an appropriately sized window.  The name of the
%  button that is pressed is returned in ButtonName.  The Title of
%  the figure may be specified by adding a second string argument:
%
%    ButtonName = questdlg(Question, Title)
%
%  Question will be interpreted as a normal string.
%
%  QUESTDLG uses UIWAIT to suspend execution until the user responds.
%
%  The default set of buttons names for QUESTDLG are 'Yes','No' and 'Cancel'.
%  The default answer for the above calling syntax is 'Yes'.
%  This can be changed by adding a third argument which specifies the
%  default Button:
%
%    ButtonName = questdlg(Question, Title, 'No')
%
%  Up to 3 custom button names may be specified by entering
%  the button string name(s) as additional arguments to the function
%  call.  If custom button names are entered, the default button
%  must be specified by adding an extra argument, DEFAULT, and
%  setting DEFAULT to the same string name as the button you want
%  to use as the default button:
%
%    ButtonName = questdlg(Question, Title, Btn1, Btn2, DEFAULT);
%
%  where DEFAULT is set to Btn1.  This makes Btn1 the default answer.
%  If the DEFAULT string does not match any of the button string names,
%  a warning message is displayed.
%
%  To use TeX interpretation for the Question string, a data
%  structure must be used for the last argument, i.e.
%
%    ButtonName = questdlg(Question, Title, Btn1, Btn2, OPTIONS);
%
%  The OPTIONS structure must include the fields Default and Interpreter.
%  Interpreter may be 'none' or 'tex' and Default is the default button
%  name to be used.
%
%  If the dialog is closed without a valid selection, the return value
%  is empty.
%
%  Example:
%
%  ButtonName = questdlg('What is your favorite color?', ...
%                        'Color Question', ...
%                        'Red', 'Green', 'Blue', 'Green');
%  switch ButtonName,
%    case 'Red',
%     disp('Your favorite color is Red');
%    case 'Blue',
%     disp('Your favorite color is Blue.')
%     case 'Green',
%      disp('Your favorite color is Green.');
%  end % switch
%
%  See also DIALOG, ERRORDLG, HELPDLG, INPUTDLG, LISTDLG,
%    MSGBOX, WARNDLG, FIGURE, TEXTWRAP, UIWAIT, UIRESUME.


%  Copyright 1984-2011 The MathWorks, Inc.



if nargin<1
    error('Too few arguments for BUTTONDLG');           % <-- Modified
end

Interpreter='none';
Question = dialogCellstrHelper(Question);
needsLookup = false;
IconString ='quest';                                    % <-- Added
IconData   =[];                                         % <-- Added
IconCMap   =[];                                         % <-- Added

%%%%%%%%%%%%%%%%%%%%%
%%% General Info. %%%
%%%%%%%%%%%%%%%%%%%%%
Black      =[0       0        0      ]/255;
% LightGray  =[192     192      192    ]/255;
% LightGray2 =[160     160      164    ]/255;
% MediumGray =[128     128      128    ]/255;
% White      =[255     255      255    ]/255;

%%%%%%%%%%%%%%%%%%%%
%%% Nargin Check %%%
%%%%%%%%%%%%%%%%%%%%
if nargout > 1
    error('Wrong number of output arguments for BUTTONDLG'); % <-- Modified
end
if nargin == 1
    Title = ' ';
end
if nargin <= 2,
    Default = 'Yes'; 
    needsLookup = true;
end
if nargin == 3,
    Default = Btn1;
end
if nargin <= 3,
    Btn1 = 'Yes'; 
    Btn2 = 'No';
    Btn3 = 'Cancel';
    NumButtons = 3;
    needsLookup = true;
end
if nargin == 4,
    Default=Btn2;
    Btn2 = [];
    Btn3 = [];
    NumButtons = 1;
end
if nargin == 5
    Default = Btn3;
    Btn3 = [];
    NumButtons = 2;
end
if nargin == 6
    NumButtons = 3;
end
if nargin > 6
    error(message('MATLAB:questdlg:TooManyInputs'));
    NumButtons = 3;
end

if isstruct(Default),
    if isfield(Default,'Interpreter')                   % <-- Added
        Interpreter = Default.Interpreter;
    end                                                 % <-- Added
    if isfield(Default,'IconString')                    % <-- Added
        IconString = Default.IconString;                % <-- Added
    end                                                 % <-- Added
    if isfield(Default,'IconData')                      % <-- Added
        IconData = Default.IconData;                    % <-- Added
    end                                                 % <-- Added
    if isfield(Default,'IconCMap')                      % <-- Added
        IconCMap = Default.IconCMap;                    % <-- Added
    end                                                 % <-- Added
    if isfield(Default,'Default')                       % <-- Added
        Default = Default.Default;
    else                                                % <-- Added
        Default = '';                                   % <-- Added
    end                                                 % <-- Added
end

% Added Begin -->
IconString=lower(IconString);
switch(IconString)
    case {'custom'}
        % check for icon data
        if isempty(IconData)
            error(message('MATLAB:msgbox:icondata'))
        end
        if ~isnumeric(IconData)
            error(message('MATLAB:msgbox:IncorrectIconDataType'))
        end
        if ~isnumeric(IconCMap)
            error(message('MATLAB:msgbox:IncorrectIconColormap'))
        end
    case {'none','help','warn','error','quest','bulb'}
        % icon String OK
    otherwise
        warning('Invalid string for Icon in BUTTONDLG.');
        IconString='none';
end                                                     % <-- Added End

%%%%%%%%%%%%%%%%%%%%%%%
%%% Create QuestFig %%%
%%%%%%%%%%%%%%%%%%%%%%%
FigPos    = get(0,'DefaultFigurePosition');
FigPos(3) = 267;
FigPos(4) =  70;
FigPos    = getnicedialoglocation(FigPos, get(0,'DefaultFigureUnits'));

QuestFig=dialog(                                    ...
    'Visible'         ,'off'                      , ...
    'Name'            ,Title                      , ...
    'Pointer'         ,'arrow'                    , ...
    'Position'        ,FigPos                     , ...
    'KeyPressFcn'     ,@doFigureKeyPress          , ...
    'IntegerHandle'   ,'off'                      , ...
    'WindowStyle'     ,'normal'                   , ...
    'HandleVisibility','callback'                 , ...
    'CloseRequestFcn' ,@doDelete                  , ...
    'Tag'             ,Title                        ...
    );

%%%%%%%%%%%%%%%%%%%%%
%%% Set Positions %%%
%%%%%%%%%%%%%%%%%%%%%
DefOffset  =10;

IconWidth  =50;                                         % <-- Modified
IconHeight =50;                                         % <-- Modified
IconXOffset=DefOffset;
IconYOffset=FigPos(4)-DefOffset-IconHeight;  %#ok
if isempty(IconCMap)                                    % <-- Added
    IconCMap=[Black;get(QuestFig,'Color')];  %#ok
end                                                     % <-- Added

DefBtnWidth =56;
BtnHeight   =22;

BtnYOffset=DefOffset;

BtnWidth=DefBtnWidth;

ExtControl=uicontrol(QuestFig   , ...
    'Style'    ,'pushbutton', ...
    'String'   ,' '          ...
    );

btnMargin=1.4;
set(ExtControl,'String',Btn1);
BtnExtent=get(ExtControl,'Extent');
BtnWidth=max(BtnWidth,BtnExtent(3)+8);
if NumButtons > 1
    set(ExtControl,'String',Btn2);
    BtnExtent=get(ExtControl,'Extent');
    BtnWidth=max(BtnWidth,BtnExtent(3)+8);
    if NumButtons > 2
        set(ExtControl,'String',Btn3);
        BtnExtent=get(ExtControl,'Extent');
        BtnWidth=max(BtnWidth,BtnExtent(3)*btnMargin);
    end
end
BtnHeight = max(BtnHeight,BtnExtent(4)*btnMargin);

delete(ExtControl);

MsgTxtXOffset=IconXOffset+IconWidth;

FigPos(3)=max(FigPos(3),MsgTxtXOffset+NumButtons*(BtnWidth+2*DefOffset));
set(QuestFig,'Position',FigPos);

BtnXOffset=zeros(NumButtons,1);

if NumButtons==1,
    BtnXOffset=(FigPos(3)-BtnWidth)/2;
elseif NumButtons==2,
    BtnXOffset=[MsgTxtXOffset
        FigPos(3)-DefOffset-BtnWidth];
elseif NumButtons==3,
    BtnXOffset=[MsgTxtXOffset
        0
        FigPos(3)-DefOffset-BtnWidth];
    BtnXOffset(2)=(BtnXOffset(1)+BtnXOffset(3))/2;
end

MsgTxtYOffset=DefOffset+BtnYOffset+BtnHeight;
% Calculate current msg text width and height. If negative,
% clamp it to 1 since its going to be recalculated/corrected later
% based on the actual msg string
MsgTxtWidth=max(1, FigPos(3)-DefOffset-MsgTxtXOffset-IconWidth);
MsgTxtHeight=max(1, FigPos(4)-DefOffset-MsgTxtYOffset);

MsgTxtForeClr=Black;
MsgTxtBackClr=get(QuestFig,'Color');

CBString='uiresume(gcbf)';
DefaultValid = false;
DefaultWasPressed = false;
BtnHandle = cell(NumButtons, 1);
DefaultButton = 0;

% Check to see if the Default string passed does match one of the
% strings on the buttons in the dialog. If not, throw a warning.
for i = 1:NumButtons
    switch i
        case 1
            ButtonString=Btn1;
            ButtonTag='Btn1';
            if strcmp(ButtonString, Default)
                DefaultValid = true;
                DefaultButton = 1;
            end
            
        case 2
            ButtonString=Btn2;
            ButtonTag='Btn2';
            if strcmp(ButtonString, Default)
                DefaultValid = true;
                DefaultButton = 2;
            end
        case 3
            ButtonString=Btn3;
            ButtonTag='Btn3';
            if strcmp(ButtonString, Default)
                DefaultValid = true;
                DefaultButton = 3;
            end
    end
    
    if (needsLookup)
        buttonDisplayString = getString(message(['MATLAB:uistring:popupdialogs:' ButtonString]));
    else
        buttonDisplayString = ButtonString;
    end
    
    BtnHandle{i}=uicontrol(QuestFig            , ...
        'Style'              ,'pushbutton', ...
        'Position'           ,[ BtnXOffset(1) BtnYOffset BtnWidth BtnHeight ]           , ...
        'KeyPressFcn'        ,@doControlKeyPress , ...
        'Callback'           ,CBString    , ...
        'String'             ,buttonDisplayString, ...
        'HorizontalAlignment','center'    , ...
        'Tag'                ,ButtonTag     ...
        );
    
    setappdata(BtnHandle{i},'QuestDlgReturnName',ButtonString);   
end

if ~DefaultValid
    warnstate = warning('backtrace','off');
    warning(message('MATLAB:questdlg:StringMismatch'));
    warning(warnstate);
end

MsgHandle=uicontrol(QuestFig            , ...
    'Style'              ,'text'         , ...
    'Position'           ,[MsgTxtXOffset MsgTxtYOffset 0.95*MsgTxtWidth MsgTxtHeight ]              , ...
    'String'             ,{' '}          , ...
    'Tag'                ,'Question'     , ...
    'HorizontalAlignment','left'         , ...
    'FontWeight'         ,'bold'         , ...
    'BackgroundColor'    ,MsgTxtBackClr  , ...
    'ForegroundColor'    ,MsgTxtForeClr    ...
    );

[WrapString,NewMsgTxtPos]=textwrap(MsgHandle,Question,75);

% NumLines=size(WrapString,1);

AxesHandle=axes('Parent',QuestFig,'Position',[0 0 1 1],'Visible','off');

texthandle=text( ...
    'Parent'              ,AxesHandle                      , ...
    'Units'               ,'pixels'                        , ...
    'Color'               ,get(BtnHandle{1},'ForegroundColor')   , ...
    'HorizontalAlignment' ,'left'                          , ...
    'FontName'            ,get(BtnHandle{1},'FontName')    , ...
    'FontSize'            ,get(BtnHandle{1},'FontSize')    , ...
    'VerticalAlignment'   ,'bottom'                        , ...
    'String'              ,WrapString                      , ...
    'Interpreter'         ,Interpreter                     , ...
    'Tag'                 ,'Question'                        ...
    );

textExtent = get(texthandle, 'Extent');

% (g357851)textExtent and extent from uicontrol are not the same. For window, extent from uicontrol is larger
%than textExtent. But on Mac, it is reverse. Pick the max value.
MsgTxtWidth=max([MsgTxtWidth NewMsgTxtPos(3)+2 textExtent(3)]);
MsgTxtHeight=max([MsgTxtHeight NewMsgTxtPos(4)+2 textExtent(4)]);

MsgTxtXOffset=IconXOffset+IconWidth+DefOffset;
FigPos(3)=max(NumButtons*(BtnWidth+DefOffset)+DefOffset, ...
    MsgTxtXOffset+MsgTxtWidth+DefOffset);


% Center Vertically around icon
if IconHeight>MsgTxtHeight,
    IconYOffset=BtnYOffset+BtnHeight+DefOffset;
    MsgTxtYOffset=IconYOffset+(IconHeight-MsgTxtHeight)/2;
    FigPos(4)=IconYOffset+IconHeight+DefOffset;
    % center around text
else
    MsgTxtYOffset=BtnYOffset+BtnHeight+DefOffset;
    IconYOffset=MsgTxtYOffset+(MsgTxtHeight-IconHeight)/2;
    FigPos(4)=MsgTxtYOffset+MsgTxtHeight+DefOffset;
end

if NumButtons==1,
    BtnXOffset=(FigPos(3)-BtnWidth)/2;
elseif NumButtons==2,
    BtnXOffset=[(FigPos(3)-DefOffset)/2-BtnWidth
        (FigPos(3)+DefOffset)/2
        ];
    
elseif NumButtons==3,
    BtnXOffset(2)=(FigPos(3)-BtnWidth)/2;
    BtnXOffset=[BtnXOffset(2)-DefOffset-BtnWidth
        BtnXOffset(2)
        BtnXOffset(2)+BtnWidth+DefOffset
        ];
end

set(QuestFig ,'Position',getnicedialoglocation(FigPos, get(QuestFig,'Units')));
assert(iscell(BtnHandle));
BtnPos=cellfun(@(bh)get(bh,'Position'), BtnHandle, 'UniformOutput', false);
BtnPos=cat(1,BtnPos{:});
BtnPos(:,1)=BtnXOffset;
BtnPos=num2cell(BtnPos,2);

assert(iscell(BtnPos));
cellfun(@(bh,pos)set(bh, 'Position', pos), BtnHandle, BtnPos, 'UniformOutput', false);

if DefaultValid
    setdefaultbutton(QuestFig, BtnHandle{DefaultButton});
end

delete(MsgHandle);


set(texthandle, 'Position',[MsgTxtXOffset MsgTxtYOffset 0]);

% create icon
if ~strcmp(IconString,'none')                               % <-- Added
    IconAxes=axes(                                      ...
        'Parent'      ,QuestFig              , ...
        'Units'       ,'Pixels'              , ...
        'Position'    ,[IconXOffset IconYOffset IconWidth IconHeight], ...
        'NextPlot'    ,'replace'             , ...
        'Tag'         ,'IconAxes'              ...
        );

    set(QuestFig ,'NextPlot','add');
    % Added Begin -->
    if ~strcmp(IconString,'custom')
        a = load('dialogicons.mat');
        switch(IconString)
            case'warn'
                IconData=a.warnIconData;
                a.warnIconMap(256,:)=get(QuestFig,'Color');
                IconCMap=a.warnIconMap;

            case 'help'
                IconData=a.helpIconData;
                a.helpIconMap(256,:)=get(QuestFig,'Color');
                IconCMap=a.helpIconMap;

            case 'error'
                IconData=a.errorIconData;
                a.errorIconMap(146,:)=get(QuestFig,'Color');
                IconCMap=a.errorIconMap;
                
            case 'quest'
                IconData=a.questIconData;
                a.questIconMap(256,:)=get(QuestFig,'Color');
                IconCMap=a.questIconMap;
                
             case 'bulb'
                IconData=a.lightbulbIconData;
                a.lightbulbIconMap(6,:)=get(QuestFig,'Color');
                IconCMap=a.lightbulbIconMap;
        end
    end

    % place the icon - if this fails, rethrow the error
    % after deleting the figure
    try                                                 % <-- Added End
        Img=image('CData',IconData,'Parent',IconAxes);
        set(QuestFig, 'Colormap', IconCMap);
    catch ex % Added Begin -->
        delete(QuestFig);
        rethrow(ex);
    end
    
    if ~isempty(get(Img,'XData')) && ~isempty(get(Img,'YData')) % <-- Added End
        set(IconAxes          , ...
            'XLim'            ,get(Img,'XData')+[-0.5 0.5], ... % <-- Modified
            'YLim'            ,get(Img,'YData')+[-0.5 0.5]  ... % <-- Modified
            );
    end

    set(IconAxes          , ...
        'Visible'         ,'off'       , ...
        'YDir'            ,'reverse'     ...
        );

end % if ~strcmp                                        % <-- Added

% make sure we are on screen
movegui(QuestFig)


set(QuestFig ,'WindowStyle','modal','Visible','on');
drawnow;

if DefaultButton ~= 0
    uicontrol(BtnHandle{DefaultButton});
end

if ishghandle(QuestFig)
    % Go into uiwait if the figure handle is still valid.
    % This is mostly the case during regular use.
    uiwait(QuestFig);
end

% Check handle validity again since we may be out of uiwait because the
% figure was deleted.
if ishghandle(QuestFig)
    if DefaultWasPressed
        ButtonName=Default;
    else
        ButtonName = getappdata(get(QuestFig,'CurrentObject'),'QuestDlgReturnName');
    end
    doDelete;
else
    ButtonName='';
end

    function doFigureKeyPress(obj, evd)  %#ok
    switch(evd.Key)
        case {'return','space'}
            if DefaultValid
                DefaultWasPressed = true;
                uiresume(gcbf);
            end
        case 'escape'
            doDelete
    end
    end

    function doControlKeyPress(obj, evd)  %#ok
    switch(evd.Key)
        case {'return'}
            if DefaultValid
                DefaultWasPressed = true;
                uiresume(gcbf);
            end
        case 'escape'
            doDelete
    end
    end

    function doDelete(varargin)
    delete(QuestFig);
    end

    % all functions beginning here are not changed, but needs to be
    % importet, no idea why
    function outStr = dialogCellstrHelper (inputStr)
    %   Copyright 2010-2011 The MathWorks, Inc.

    % Helper used by MSGBOX, ERRORDLG, WARNDLG, QUESTDLG to parse the input
    % string vector, matrix or cell array or strings.
    % This works similar to the CELLSTR function but does not use deblank, like
    % cellstr, to eliminate any trailing white spaces.

    % Validate input string type. 
    validateattributes(inputStr, {'char','cell'}, {'2d'},mfilename);

    % Convert to cell array of strings without eliminating any user input. 
    if ~iscell(inputStr)
        inputCell = {};
        for siz = 1:size(inputStr,1)
            inputCell{siz} =inputStr(siz,:); %#ok<AGROW>
        end
        outStr = inputCell;
    else
        outStr = inputStr;
    end
    end % end function dialogCellstrHelpber

    function figure_size = getnicedialoglocation(figure_size, figure_units)
    % adjust the specified figure position to fig nicely over GCBF
    % or into the upper 3rd of the screen

    %  Copyright 1999-2010 The MathWorks, Inc.
    %  $Revision: 1.1.6.5 $

    parentHandle = gcbf;
    convertData.destinationUnits = figure_units;
    if ~isempty(parentHandle)
        % If there is a parent figure
        convertData.hFig = parentHandle;
        convertData.size = get(parentHandle,'Position');
        convertData.sourceUnits = get(parentHandle,'Units');  
        c = []; 
    else
        % If there is no parent figure, use the root's data
        % and create a invisible figure as parent
        convertData.hFig = figure('visible','off');
        convertData.size = get(0,'ScreenSize');
        convertData.sourceUnits = get(0,'Units');
        c = onCleanup(@() close(convertData.hFig));
    end

    % Get the size of the dialog parent in the dialog units
    container_size = hgconvertunits(convertData.hFig, convertData.size ,...
        convertData.sourceUnits, convertData.destinationUnits, get(convertData.hFig,'Parent'));

    delete(c);

    figure_size(1) = container_size(1)  + 1/2*(container_size(3) - figure_size(3));
    figure_size(2) = container_size(2)  + 2/3*(container_size(4) - figure_size(4));
    end % end function getnicedialoglocation

    function setdefaultbutton(figHandle, btnHandle)
    % WARNING: This feature is not supported in MATLAB and the API and
    % functionality may change in a future release.

    %SETDEFAULTBUTTON Set default button for a figure.
    %  SETDEFAULTBUTTON(BTNHANDLE) sets the button passed in to be the default button
    %  (the button and callback used when the user hits "enter" or "return"
    %  when in a dialog box.
    %
    %  This function is used by inputdlg.m, msgbox.m, questdlg.m and
    %  uigetpref.m.
    %
    %  Example:
    %
    %  f = figure;
    %  b1 = uicontrol('style', 'pushbutton', 'string', 'first', ...
    %       'position', [100 100 50 20]);
    %  b2 = uicontrol('style', 'pushbutton', 'string', 'second', ...
    %       'position', [200 100 50 20]);
    %  b3 = uicontrol('style', 'pushbutton', 'string', 'third', ...
    %       'position', [300 100 50 20]);
    %  setdefaultbutton(b2);
    %

    %  Copyright 2005-2007 The MathWorks, Inc.

    %--------------------------------------- NOTE ------------------------------------------
    % This file was copied into matlab/toolbox/local/private.
    % These two files should be kept in sync - when editing please make sure
    % that *both* files are modified.

    % Nargin Check
    error(nargchk(1, 2, nargin))

    if (usejava('awt') == 1)
        % We are running with Java Figures
        useJavaDefaultButton(figHandle, btnHandle)
    else
        % We are running with Native Figures
        useHGDefaultButton(figHandle, btnHandle);
    end

        function useJavaDefaultButton(figH, btnH)
            % Get a UDD handle for the figure.
            fh = handle(figH);
            % Call the setDefaultButton method on the figure handle
            fh.setDefaultButton(btnH);
        end

        function useHGDefaultButton(figHandle, btnHandle)
            % First get the position of the button.
            btnPos = getpixelposition(btnHandle);

            % Next calculate offsets.
            leftOffset   = btnPos(1) - 1;
            bottomOffset = btnPos(2) - 2;
            widthOffset  = btnPos(3) + 3;
            heightOffset = btnPos(4) + 3;

            % Create the default button look with a uipanel.
            % Use black border color even on Mac or Windows-XP (XP scheme) since
            % this is in natve figures which uses the Win2K style buttons on Windows
            % and Motif buttons on the Mac.
            h1 = uipanel(get(btnHandle, 'Parent'), 'HighlightColor', 'black', ...
                'BorderType', 'etchedout', 'units', 'pixels', ...
                'Position', [leftOffset bottomOffset widthOffset heightOffset]);

            % Make sure it is stacked on the bottom.
            uistack(h1, 'bottom');
        end % end function useJavaDefaultButton
    end % end function setdefaultbutton

end % end function buttondlg