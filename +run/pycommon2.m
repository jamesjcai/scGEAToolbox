function [status] = pycommon2(x, wrkpth, prgname)

fw = gui.myWaitbar([], [], [], sprintf('Running %s...', ...
upper(strrep(prgname, '_', '\_'))));
cmdlinestr = sprintf('"%s" "%s%sscript.py"', ...
x.Executable, wrkpth, filesep);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');
if pkg.i_isvalid(fw)
    gui.myWaitbar([], fw, [], ...
        sprintf('Running %s is complete', ...
        upper(strrep(prgname, '_', '\_'))));
end
end
