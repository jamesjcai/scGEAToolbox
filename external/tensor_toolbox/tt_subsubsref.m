function a = tt_subsubsref(obj, s)
%TT_SUBSUBSREF Helper function for tensor toolbox subsref.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if isscalar(s)
    a = obj;
else
    a = subsref(obj, s(2:end));
end
