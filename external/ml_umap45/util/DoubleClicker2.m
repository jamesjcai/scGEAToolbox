classdef DoubleClicker2 < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties
        axes=0;
        multiSelectWhenButtonDown=false;
        doubleClickDelay=.33;
        
    end
    properties(SetAccess=private, GetAccess=private)
        chk=[];
        cp=[];
        figure=0;
        fncSingleUp=0;
        fncFirstDown=0;
        fncSingleDown=0;
        fncDoubleUp=0;
        fncDoubleDown=0;
        state=0;
        ups=0;
        downs=0;
        tmr;
    end
    
    methods
        function this=DoubleClicker2(figure, axes)
            this.figure=figure;
            this.axes=axes;
            this.startListening();
            this.tmr=tic;
        end
        
        function setFnc(this, fncSingleDown,fncSingleUp, ...
                fncDoubleDown, fncDoubleUp, fncFirstDown)
            if nargin<2
                fncSingleUp=0;
            end
            if nargin<3
                fncSingleDown=0;
            end
            if nargin<4
                fncDoubleUp=0;
            end
            if nargin<5
                fncDoubleDown=0;
            end
            this.fncDoubleUp=fncDoubleUp;
            this.fncDoubleDown=fncDoubleDown;
            this.fncSingleUp=fncSingleUp;
            this.fncSingleDown=fncSingleDown;
            this.fncFirstDown=fncFirstDown;
        end
        
        function down(this, obj,evt)
            this.multiSelectWhenButtonDown=DoubleClicker2.isMultiSelect;
            if ~isempty(this.axes)
                this.cp=get(this.axes, 'CurrentPoint');
            end
            if ~ishandle(obj)
                return;
            end
            this.downs=this.downs+1;
                fprintf('nice');
                secs=toc(this.tmr);
                disp(secs)
            if secs>this.doubleClickDelay%isempty(this.chk)
                if ~isnumeric(this.fncFirstDown)
                    if feval(this.fncFirstDown, obj, this.cp, guidata(obj))
                        return;
                    end
                end
                this.chk=1;
                this.state=0;
                if this.chk==1
                    this.state=1;
                    this.chk=[];
                    if ~isnumeric(this.fncSingleDown)
                        feval(this.fncSingleDown, obj, this.cp, guidata(obj));
                    end
                    if this.ups==this.downs || this.downs<2
                        %fprintf(1,'\nButton DOWN single-click WITH button UP.\n\n');
                        this.up(obj,evt);
                    else
                        %fprintf(1,'\nButton DOWN single-click with NO button UP.\n\n');
                    end
                end
            else
                this.state=2;
                this.chk=[];
                if ~isnumeric(this.fncDoubleDown)
                    feval(this.fncDoubleDown, obj, this.cp, guidata(obj));
                end
                %fprintf(1,'\nButton DOWN a double-click.\n\n');
            end
            this.tmr=tic;
        end
        
        function up(this, obj,evt)
            if ~ishandle(obj)
                return;
            end
            
            this.ups=this.downs;
            if this.state==0
                fprintf(1,'\nButton UP waiting for double or single click.\n\n');
            elseif this.state==1
                %fprintf(1,'\nButton UP a single-click.\n\n');
                if ~isnumeric(this.fncSingleUp)
                    feval(this.fncSingleUp, obj, this.cp, guidata(obj));
                end
                
            elseif this.state==2
                %fprintf(1,'\nButton UP a double-click.\n\n');
                if ~isnumeric(this.fncDoubleUp)
                    feval(this.fncDoubleUp, obj, this.cp, guidata(obj));
                end
                
            end
        end
        
        function stopListening(this)
            set(this.figure,'WindowButtonDownFcn', '');
            set(this.figure, 'WindowButtonUpFcn', '');
        end
        
        function startListening(this)            
            set(this.figure,'WindowButtonDownFcn', @(obj,evt)down(this,obj,evt));
            set(this.figure, 'WindowButtonUpFcn', @(obj,evt)up(this,obj,evt));
        end
        
        

    end
    
    methods(Static)
        function ok=isMultiSelect
            if ismac
                ok=edu.stanford.facs.swing.KeyState.IsMetaPressed;
            else
                ok=edu.stanford.facs.swing.KeyState.IsCtrlPressed;
            end
        end

        function test()
            treeTest;
            DoubleClicker2(gcf);
        end
    end
end
