classdef SuhAbstractClass< handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause

    methods(Static)
        function AssertIsA(this, className)
            assert(isa(this, className),...
                'set arg must be instance of %s', className);

        end
    end
    
    properties(SetAccess=private)
        sac_createdWhen;
        sac_id;
        sac_updatedWhen;
    end
    
    methods(Sealed)
        function this=SuhAbstractClass()
           this.sac_createdWhen=now;
        end
        
        function set_sac_id(this, id)
            this.sac_id=id;
        end
        
        function warnNotImplemented(this)
            funCallStack = dbstack;
            if length(funCallStack)>1
                methodName = funCallStack(2).name;
            else
                methodName='commandConsole';
            end
            warning('%s has not implemented function %s()', ...
                class(this), methodName);
        end
    end
    
    methods
        function noteUpdated(this)
            this.sac_updatedWhen=now;
        end
        
        function ok=describeProgress(~) %this
            ok=true;
        end
                    
    end
end