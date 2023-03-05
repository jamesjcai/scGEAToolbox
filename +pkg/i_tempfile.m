function [b,a]=i_tempfile(subdirname)
if nargin<1
    subdirname=[];
end
if ~isempty(subdirname)
    a=fullfile(tempdir,subdirname);
    if ~exist(a,"dir")
        mkdir(a);
    end
else
    a=tempdir;
end

tstr=matlab.lang.makeValidName(string(datetime));
b=fullfile(a,strcat(tstr,".txt"));
end