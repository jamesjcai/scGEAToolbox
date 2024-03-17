classdef SuhWebUpdate < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
     
    properties(SetAccess=private)
        map;
        appName;
        exeFile;
        urlPing;
        urlMac;
        urlPc;
        urlVersionFile;
        urlSourceCode;
        urlHome;
        localVersion;
        serverVersion;
        v;
    end
    
    methods
        function this=SuhWebUpdate(app)
            [this.v, this.localVersion, this.appName, this.exeFile,...
                versionFile, remoteFolder, remoteSourceCode]...
                =app.getAppDetails;
            app.rememberUpdateCheck;
            url=WebDownload.ResolveUrl(versionFile, remoteFolder);
            host=WebDownload.UrlParts(url);
            if isempty(host)
                this.urlPing='**no host found**';
            else
                this.urlPing=host;
            end
            remotePath=fileparts(url);
            J=edu.stanford.facs.swing.WebDownload;
            this.urlMac=J.GoogleDirectDownloadLink(...
                fullfile(remotePath, this.exeFile));
            this.urlPc=this.urlMac;
            this.urlVersionFile=WebDownload.FullUrl(remotePath, versionFile);
            if ~isempty(remoteSourceCode)
                try
                    %uploaded umapAndEpp.zip to Google Cloud on March 29,
                    %2022
                    this.urlSourceCode=fullfile(remotePath(...
                        1:String.IndexOf(remotePath, remoteFolder)-1), ...
                        remoteSourceCode);
                catch
                end
            end
            this.urlHome='http://www.cytogenie.org';
        end
        
        function ok=isReachable(this)
            ok=edu.stanford.facs.swing.CpuInfo.isReachable(...
                this.urlPing, 80, 5500);
        end
        
        function [yesNewVersion, connectedOk, serverVersion, ...
                localVersion]=hasNew(this)
            connectedOk=this.isReachable;
            if connectedOk
                try
                    serverVersion=str2double(WebDownload.ReadText(this.urlVersionFile));
                    yesNewVersion=serverVersion>this.localVersion;
                catch
                    serverVersion=nan;
                    yesNewVersion=false;
                end
            else
                yesNewVersion=false;
                serverVersion=this.localVersion;
            end
            localVersion=this.localVersion;
            this.serverVersion=serverVersion;
        end
        
        function tip=tipText(this, showMsg)
            if ismac
                finderName='Mac finder';
                [~, fileName, ext]=fileparts(this.exeFile);
                fileNameCopy=[fileName '-2' ext];
                li3=['<li>Open <u>' finderName '</u> on the folder & double click '...
                    ' <b>' fileName ext '</b> to <br>&nbsp;&nbsp;unzip <b>' ...
                    fileName '.app</b>'];
                li4=['<li>Double click  <b>' fileName ...
                    '.app</b> to run ' this.appName];
            else
                finderName='File explorer';
                [~,fileName, ext]=fileparts(this.urlPc);
                fileNameCopy=[fileName '- Copy' ext];
                li3=['<li>Open <u>' finderName '</u> on the folder & double click '...
                    ' <b>' fileName ext '</b> to run ' this.appName];
                li4='';
            end
            tip=['After downloading <ol>'...
                '<li>Move <b>' fileName ext '</b> to a folder '...
                    'from where you want to run AutoGate'...
                '<li>If file of same name pre-exists, this '...
                    'renames it as <b>' fileNameCopy '</b>'...
                li3 li4 '.</ol><hr>'];

            if nargin>1 && showMsg
                msg(Html.Wrap(tip), 40);
            end
        end
        
        function reportToUser(this, connectedOk)
            if connectedOk
                msg(Html.WrapHr(...
                    ['<b>' this.appName '</b> has no updates...<br><br>'...
                    'You are <b>up to date</b>!']), 7, 'south');
            else
                msg(Html.Wrap(['Could not connect to the server <b>' ...
                    this.urlPing '</b> <br>Re-check for updates later.']));
            end
        end
        
        function ask(this, userIsAsking)
            if userIsAsking
                pauseSecs=0;
            else
                pauseSecs=7;
            end
            questDlg(struct('msg', Html.WrapHr(...
                ['<b>The newer version ' num2str(this.serverVersion) ...
                ' is available</b>...<br><i>(Your current version is '...
                num2str(this.localVersion) ')</i><br><br>'...
                Html.WrapSmallTags(['(Click <u>Download</u> button ' ...
                ' to download version <b>' num2str(this.serverVersion) ...
                '</b><br>... the local version is <b>' ...
                num2str(this.localVersion) '</b>)'])]), ...
                'checkFnc', @(jd, answer)download(jd, answer), ...
                'modal', false, 'pauseSecs', pauseSecs), ...
                ['Download new ' this.appName '?'], ...
                'Download', 'Ok', 'Download');
            
            
            function ok=download(jd, answer)
                ok=false;
                if isequal('Download', answer)
                    isGoogleDriveDir=startsWith(this.urlVersionFile, ...
                            'https://drive.google.com');

                    if ~isdeployed && strcmpi(this.appName, 'FlowJoBridge')
                        if exist(fullfile(fileparts(fileparts( ...
                                fileparts(mfilename('fullpath')))), ...
                                '.git'), 'dir')
                            if ~askYesOrNo(struct('javaWindow', jd,...
                                    'msg', Html.WrapHr(['All of this ' ...
                                    'source code is in <i>GitHub</i>... ' ...
                                    '<br>Download <i>anyway</i>?'])), 'Confirm...', 'north++')
                                return;
                            end
                        end
                        to=File.Downloads('UmapAndFlowJoBridge.zip');
                        [cancelled, bad]=WebDownload.Get({this.urlSourceCode}, to, false);
                        if cancelled
                            return;
                        end
                        if ~bad
                            SuhWebUpdate.TryNewSourceCode(to, jd);
                        else
                            msgError(struct('javaWindow', jd, 'msg', ...
                                Html.Wrap(['Could not ' ...
                                'download' Html.FileTree(url)])));
                        end
                    else
                        if isGoogleDriveDir
                            Html.BrowseString(tipText(this, false));
                        else
                            web(this.urlHome, '-browser');
                        end
                        pause(1);
                        if ismac
                            web(this.urlMac, '-browser');
                        else
                            web(this.urlPc, '-browser');
                        end
                        if ~isempty(this.urlSourceCode)
                            if askYesOrNo('Download source too?')
                                web(this.urlSourceCode, '-browser');
                            end
                        end
                        MatBasics.RunLater(@(h,e)tipText(this, true), 3);
                    end
                end
                ok=true;
                jd.dispose;                
            end
        end
    end
    
    methods(Static)
        function [madeConnection, app]=Check(app, userIsAsking, waitSecs, ...
                helpNow, helpBefore, propTimeUpdate, hoursBetweenChecks)
            if nargin<7
                hoursBetweenChecks=48;
                if nargin<6
                    propTimeUpdate=[];
                    if nargin<5
                        helpBefore=[];
                        if nargin<4
                            helpNow=[];
                            if nargin<3
                                waitSecs=.25;
                                if nargin<2
                                    userIsAsking=false;
                                end
                            end
                        end
                    end
                end
            end
            madeConnection=false;
            if userIsAsking
                pu=PopUp('Checking cloud for updates');
                Gui.Modem;
                if ~isempty(helpNow)
                    BasicMap.Global.setHelp(helpNow);
                end
                check;
                if ~isempty(helpBefore)
                    BasicMap.Global.setHelp(helpBefore);
                end
                pu.close;
            else
                if ~isempty(propTimeUpdate)
                    t1=BasicMap.Global.getNumeric(propTimeUpdate, nan);
                    if isnan(t1) 
                        since=hoursBetweenChecks*2;
                    else
                        since=minutes(diff(datetime(...
                            [datevec(t1);datevec(now)])))/60;
                    end
                    if since<hoursBetweenChecks
                        return;
                    end
                end
                pu=PopUp('Checking cloud for updates');
                MatBasics.RunLater(@(h,e)check, waitSecs);
            end
            
            function check
                wu=SuhWebUpdate(app);
                [yesNewVersion, madeConnection]=wu.hasNew;
                if madeConnection
                    if ~isempty(propTimeUpdate)
                        BasicMap.Global.setNumeric(propTimeUpdate,...
                            num2str(now));
                    end
                    if yesNewVersion
                        pu.close;
                        wu.ask(userIsAsking);
                        return;
                    end                    
                end
                pu.close;
                    
                if userIsAsking
                    wu.reportToUser(madeConnection);
                end
            end
        end

        function TryNewSourceCode(to, jd)
            [~,f,e]=fileparts(to);
            if askYesOrNo(struct('javaWindow', jd, 'msg', Html.Wrap( ...
                    ['To try out our new source code:<ul><li>Unzip <u>' ...
                    f e '</u> into a new folder.<li>Restart MATLAB.' ...
                    '<li>Change MATLAB''s current folder to the ' ...
                    'unzip folder.<li>Type the command <b>flb</b>.</ul>' ...
                    '<hr><b>Open the folder that contains <u>' f e ...
                    '</u>?'])), 'How to try our update...', 'north++', ...
                    true, 'SuhWebUpdate.OpenZip')
                File.OpenFolderWindow(to,'',false);
            end
        end

    end
end