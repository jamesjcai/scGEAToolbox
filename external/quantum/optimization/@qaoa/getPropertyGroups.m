function propgrp = getPropertyGroups(obj)
%GETPROPERTYGROUPS Return the property groups for the display
%
%   GROUPS = GETPROPERTYGROUPS(OBJ) returns the property groups for the
%   display. 
%
%   This method must be implemented as this class inherits from
%   matlab.mixin.CustomDisplay.

%   Copyright 2024 The MathWorks, Inc.

if ~isscalar(obj)
    propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
    return
end

[StrongStartTag, StrongEndTag] = iCreateStrongTags;

% Optimization properties group
optimList = ["InitialAngles", "OptimizationSolver", "OptimizationSolverOptions"];
optimHeader = getString(message("quantum:annealing:qaoa:OptimizationPropertiesHeader"));
optimTitle = StrongStartTag + string(optimHeader) + StrongEndTag;
optimGrp = matlab.mixin.util.PropertyGroup(optimList, optimTitle);

% Circuit properties group
circuitList = ["NumLayers", "NumShots"];
circuitHeader = getString(message("quantum:annealing:qaoa:CircuitPropertiesHeader"));
circuitTitle = StrongStartTag + string(circuitHeader) + StrongEndTag;
circuitGrp = matlab.mixin.util.PropertyGroup(circuitList,circuitTitle);
propgrp = [optimGrp,circuitGrp];


function [StrongStartTag, StrongEndTag] = iCreateStrongTags

if matlab.internal.display.isHot
    StrongStartTag = '<strong>';
    StrongEndTag = '</strong>';
else
    StrongStartTag = '';
    StrongEndTag = '';
end