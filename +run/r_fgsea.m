function [s] = r_fgsea(genelist, rmribo, drdistin, wkdir)

if nargin < 4, wkdir = tempdir; end
if nargin < 3, drdistin = []; end
if nargin < 2, rmribo = true; end
s=[];

isdebug = false;

oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_fgsea');
if ~isok
    error(msg);
    return;
end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.txt', 'output.txt'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

%if exist('input.txt', 'file'), delete('input.txt'); end
%if exist('output.txt', 'file'), delete('output.txt'); end

%t=readtable('input_template.txt');
%N=size(t,1);
%if length(genelist)>N
%    genelist=genelist(1:N);
%end

genelist = upper(genelist);
%a=-log(1-rand(length(genelist),1));
% a=12+randn(length(genelist),1);
% drdist=sort(a,'descend');

if ~isempty(drdistin) && length(genelist) == length(drdistin)
    drdist = drdistin;
else    
    [pw1, pw0] = cdgea;
    cd(pw0);
    dfile=fullfile(pw1,'resources','value_template_pos.txt');
    v = readmatrix(dfile);
    N = min([length(v), length(genelist)]);
    genelist = genelist(1:N);
    drdist = v(1:N);
end

sortid = (1:length(genelist))';
T = table(sortid, genelist, drdist);

if rmribo
    [gribo] = pkg.i_get_ribosomalgenes;
    i = ~ismember(T.genelist, gribo);
    T = T(i, :);
end
writetable(T, 'input.txt');
Rpath = getpref('scgeatoolbox', 'rexecutablepath');
codefullpath = fullfile(codepath,'script.R');
pkg.RunRcode(codefullpath, Rpath);

pause(1);
if exist('output.txt', 'file')
    s = readtable('output.txt', "Delimiter", ',');
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end