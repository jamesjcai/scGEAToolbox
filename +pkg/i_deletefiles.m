function i_deletefiles(a)
for k = 1:length(a)
    if exist(a{k}, 'file')
        delete(a{k});
    end
end
end