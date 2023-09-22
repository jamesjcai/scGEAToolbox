a = dir('../*export*.m');
b = struct2cell(a);

%%

names = dependencies.toolboxDependencyAnalysis(b(1, :));

pause
for k = 1:length(b(1, :))
    %    b{1,k};
    names = dependencies.toolboxDependencyAnalysis(b(1, k));
end