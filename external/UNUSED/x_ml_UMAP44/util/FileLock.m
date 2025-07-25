%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu
%   Math Lead & Secondary Developer: Connor Meehan <connor.gw.meehan@gmail.com>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
%
% based on https://stackoverflow.com/questions/3451343/atomically-creating-a-file-lock-in-matlab-file-mutex
%
classdef FileLock < handle
    properties (Access = private)
        fileName;
        fileLock=[];
        file;
        useAssociateFile;
        associateFile;
        name;
        locked=false;
    end

    methods
        %INPUT
        %fileName           To lock
        %useAssociateFile   Do NOT use Java NIO option of FileLock.m
        %                       a) Mac don't see it
        %                       b) FlowJo opens empty workspace
        %name               For display.
        function this = FileLock(fileName, useAssociateFile, name)
            if nargin<3
                name='workspace';
                if nargin<2
                    useAssociateFile=true;
                end
            end
            this.name=name;
            this.fileName=fileName;
            this.useAssociateFile=useAssociateFile;
        end

        function [missing, fileName]=lockAssociateFile(this)
            missing=false;
            if this.useAssociateFile
                fileName=File.SwitchExtension2(this.fileName, '.lock');
                while exist(fileName, 'dir')
                    fileName=[fileName '2'];
                end
                if ~exist(fileName, 'file')
                    missing=true;
                    File.WriteTextFile(fileName, 'lock');
                end
            else
                fileName='';
            end
            this.associateFile=fileName;
        end

        function ok=tryLock(this, ask)
            if nargin<2
                ask=true;
            end
            if this.useAssociateFile
                ok=this.lockAssociateFile;
                if ~ok
                    warning('%s lock file exists');
                end
            elseif isempty(this.file)
                ok=false;
                try
                    this.file = java.io.RandomAccessFile(this.fileName,'rw');
                    fileChannel = this.file.getChannel();
                    this.fileLock = fileChannel.tryLock();
                    ok=true;
                catch %ex %we know WHY ex is thrown
                    %ex.getReport
                    this.file.close;
                    this.file=[];
                end
            end
            this.locked=ok;
            if ~ok && ask
                Gui.Splat;
                ok=askYesOrNo(struct('icon', 'warning.png', 'msg',...
                    Html.WrapHr([ '<font color="red">Sorry </font>...'...
                    'another process is already<br>' ...
                    'working on this <b>' this.name '</b>!!' ...
                    '<br><br><b>Continue</b>?<br><br>' ...
                    Html.WrapSm(['(If continuing try to close other ' ...
                    'process <br>or else corruption could occur)'])])), ...
                    [String.Capitalize(this.name) ' conflict ...'], ...
                    'center', false);
            end
        end
        
        function ok = hasLock(this)
            if this.useAssociateFile
                ok=this.locked;
            elseif ~isempty(this.fileLock) && this.fileLock.isValid()
                ok = true;
            else
                ok = false;
            end
        end

        function release(this)
            if this.useAssociateFile
                f=this.associateFile;
                if ~exist(f, 'file')
                    warning('%s does not exist', f);
                else
                    delete(f);
                end
            elseif ~isempty(this.file)
                if this.hasLock
                    this.fileLock.release();
                end
                this.file.close;
                this.file=[];
            end
            this.locked=false;
        end
    end

    methods(Static)
        function MsgLock(ex)
            try
                if nargin>0 && ...
                        ~contains(ex.message, 'another process has locked')
                    Gui.MsgException(ex);
                end
            catch
                return;
            end
            msgError(Html.WrapHr(['<font color="red">Sorry</font>' ...
                '...another process is already<br>' ...
                'working on this <b>' this.name ' </b>!!']), ...
                10, 'north', [String.Capitalize(this.name) ...
                ' conflict ...']);
        end
    end
end