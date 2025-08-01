function [needupdate, v1local, v2web, im] = i_majvercheck(needed)

if nargin<1, needed = [true true true true]; end
% major version update check
needupdate = false;
v1local = [];
v2web = [];
im = [];
try    
    % fid = fopen(xfilelocal, 'r');
    % C = textscan(fid, '%s', 'delimiter', '\n');
    % fclose(fid);
    % a = C{1};    
    % x = a(contains(a, '<param.version>'));
    % a1 = strfind(x, '<param.version>');
    % a2 = strfind(x, '</param.version>');
    % v1local = extractBetween(x, a1{1}+length('<param.version>'), a2{1}-1)
    if needed(2)
        v1local = pkg.i_get_versionnum;
    end
    if needed(3)
        % xfile = 'scGEAToolbox.prj';
        % url = sprintf('https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/main/%s',xfile);
        % a = webread(url);
        % a = strsplit(a, '\n')';
        % x = a(contains(a, '<param.version>'));
        % a1 = strfind(x, '<param.version>');
        % a2 = strfind(x, '</param.version>');
        % v2web = extractBetween(x, a1{1}+length('<param.version>'), a2{1}-1);

        instURL = 'https://api.github.com/repos/jamesjcai/scGEAToolbox/releases/latest';
        instRes = webread(instURL);
        v2web = instRes.tag_name(2:end);


    end
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
    needupdate = ~isequal(v1local, v2web);
catch ME
    disp(ME.message);
end
    % if nargout > 1 && needed(2), v1local = v1local; end

    % if nargout > 2 && needed(3)
    %     if ~isempty(v2web)
    %         try
    %             v2web = v2web;
    %         catch
    %             v2web = [];
    %         end
    %     end
    % end
    if nargout > 3 && needed(4)
        try
            im = webread('https://visit-counter.vercel.app/counter.png?page=https%3A%2F%2Fgithub.com%2Fjamesjcai%2FscGEAToolbox%2F&s=15&c=ffffff&bg=00000000&no=2&ff=digi&tb=&ta=');
        catch
        end
    end
end
