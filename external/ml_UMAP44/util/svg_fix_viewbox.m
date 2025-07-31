% SVG_FIX_VIEWBOX   Make SVG file stretch with containing HTML element
%   SVG_FIX_VIEWBOX(IN_NAME,OUT_NAME) fixes the file IN_NAME so that the canvas
%   stretches to the width of the containing HTML element, when included in a
%   webpage.
%   If OUT_NAME is not given, '_fix' is appended to IN_NAME.
%
%   SVG_FIX_VIEWBOX(...,WIDTH) makes the canvas stretch along WIDTH% of the
%   containing HTML element.
%
%   Example:
%     plot(sin(1:0.1:10));
%     print -dsvg myplot.svg
%     svg_fix_viewbox myplot.svg
%     svg_fix_viewbox('myplot.svg',30,'myplot_small.svg')
%
%   Note:
%   The default figure font is written to the SVG file as 'SansSerif'. In my
%   browser (FireFox) this shows up as the standard serif font, because
%   'SansSerif' is not a known font. Correct CSS for the default sans serif font
%   is 'sans-serif'. This script will correct these fonts.

% Copyright (c)2015 by Cris Luengo
% Centre for Image Analysis, Uppsala, Sweden

function svg_fix_viewbox(in_name,varargin)

% Parameter parsing
if ~ischar(in_name)
   error('File name should be a string')
end
width = 100;
out_name = '';
for ii=1:length(varargin)
   if ischar(varargin{ii})
      if ~isempty(out_name)
         error('Unexpected input.')
      end
      out_name = varargin{ii};
   else
      if width~=100
         error('Unexpected input.')
      end
      width = varargin{ii};
      if ~isnumeric(width) || ~isscalar(width)
         error('WIDTH parameter expected to be a scalar value')
      end
   end
end


% Read SVG file
f = fopen(in_name,'rt');
if f<1
   error('File couldn''t be opened for reading');
end
s = fread(f,'*char')';
fclose(f);


% Replace 'width="xxx" height="yyy"' with 'width="100%" viewBox="0 0 xxx yyy"'
t = regexp(s,'<svg.* width="(?<width>[0-9]*)" height="(?<height>[0-9]*)"','names');
s = regexprep(s,'(?<=<svg[^\n]*) width="[0-9]*" height="[0-9]*"',sprintf(' width="%d\\%%" viewBox="0 0 %s %s"',width,t.width,t.height));


% Fix fonts -- OPTIONAL
% ('SansSerif' doesn't work on my browser, the correct CSS is 'sans-serif').
s = regexprep(s,'font-family:SansSerif;|font-family:''SansSerif'';','font-family:''sans-serif'';');
% (The document-wide default font is 'Dialog'. What is this any why?)
s = regexprep(s,'font-family:''Dialog'';','font-family:''sans-serif'';');


% Write SVG file
if isempty(out_name)
   [p,f,x] = fileparts(in_name);
   out_name = fullfile(p,[f,'_fix',x]);
end
f = fopen(out_name,'wt');
if f<1
   error('File couldn''t be opened for writing');
end
fprintf(f,'%s',s);
fclose(f);
