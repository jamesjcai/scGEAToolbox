function [needupdate,v1,v2]=i_vercheck
needupdate=false;
olddir=pwd();
cdgea;

try
%    a=textread('info.xml','%s','delimiter','\n');
fid = fopen('info.xml','r');
C = textscan(fid,'%s','delimiter','\n');
fclose(fid);
a = C{1};
% https://www.mathworks.com/matlabcentral/answers/359034-how-do-i-replace-textread-with-textscan
    x=a(contains(a,'<name>'));
    a1=strfind(x,'<name>');
    a2=strfind(x,'</name>');
    v1=extractBetween(x,a1{1}+length('<name>'),a2{1}-1);

    url='https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/master/info.xml';
    a=webread(url);
    a=strsplit(a,'\n')';
    x=a(contains(a,'<name>'));
    a1=strfind(x,'<name>');
    a2=strfind(x,'</name>');
    v2=extractBetween(x,a1{1}+length('<name>'),a2{1}-1);

    %{
    a=textread('scGEAToolbox.prj','%s');
    x=a(contains(a,'<param.version>'));
    a1=strfind(x,'<param.version>');
    a2=strfind(x,'</param.version>');
    v1=extractBetween(x,a1{1}+15,a2{1}-1);

    
    url='https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/master/scGEAToolbox.prj';
    a=webread(url);
    a=strsplit(a,'\n')';
    x=a(contains(a,'<param.version>'));
    a1=strfind(x,'<param.version>');
    a2=strfind(x,'</param.version>');
    v2=extractBetween(x,a1{1}+15,a2{1}-1);
    %}
    needupdate=~isequal(v1,v2);
catch
   
end
if nargout>1, v1=strrep(v1{1},'scGEAToolbox ',''); end
if nargout>2, v2=strrep(v2{1},'scGEAToolbox ',''); end
cd(olddir);
end
