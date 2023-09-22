classdef refwrap < handle
    %REFWRAP - a generic handle class for arbitrary MATLAB data
    %
    %  obj=refwrap(X)
    %
    %where X is any MATLAB variable will result in a handle object
    %obj with obj.data=X
    %
    %This is useful, for example, if we want to force X to be processed by
    %reference in a function call, i.e.,
    %
    %  obj=refwrap(X); clear X
    %  func(obj,...)
    %
    %would allow func() to process obj.data arbitrarily, without making a 2nd deep
    %copy of X, and so that obj.data  in the base workspace would feel the
    %changes made by func.


    properties

        data;

    end

    methods

        function obj = refwrap(dataInput)

            if nargin == 0, return; end

            obj.data = dataInput;

        end
    end

end