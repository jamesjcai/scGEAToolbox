% function [out] = x_arguments_template(X,n,b,options,axesProps)
%
% %
% web('https://www.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html')
% arguments
%     X double {mustBeNonsparse}
%     A (1,1) string
%     B (1,:) double
%     C (2,2) cell
%     n (1,1) {mustBeNumeric,mustBeReal} = 3
%     b (1,1) logical = false
%     xfinal {mustBePositive}
%     x (1,:) {mustBeNumeric,mustBeReal}
%     v (1,:) {mustBeNumeric,mustBeReal,mustBeEqualSize(v,x)}
%     method (1,:) char {mustBeMember(method,{'linear','cubic','spline'})} = 'linear'
%    options.LineStyle (1,1) string = "-"
%    options.LineWidth (1,1) {mustBeNumeric} = 1
%     axesProps.XLim (1,2) {mustBeNumeric}
%     axesProps.Color string {mustBeMember(axesProps.Color,["red","green","blue"])}
%     axesProps.Box matlab.lang.OnOffSwitchState
% end
%
% % a = myAxes("XLim",[1,10],"Color","red","Box","on");
% % myRectangle(4,5)
% % myRectangle(4,5,LineStyle=":",LineWidth=2)
% % myRectangle(4,5,LineWidth=2,LineStyle=":")
% % myRectangle(4,5,LineStyle=":")
% % myRectangle(4,5,LineWidth=2)
%
% arguments (Output)
%     out (1,:) double
% end
%
% % Custom validation function
% function mustBeEqualSize(a,b)
%     % Test for equal size
%     if ~isequal(size(a),size(b))
%         eid = 'Size:notEqual';
%         msg = 'Size of first input must equal size of second input.';
%         throwAsCaller(MException(eid,msg))
%     end
% end
%
% end
%
%
%
%
%
% https://www.mathworks.com/matlabcentral/fileexchange/59561-spider_plot?focused=68d60cc9-ef7d-474b-ba34-929e593014e4&tab=example
