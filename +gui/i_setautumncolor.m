function i_setautumncolor(c)
    a = colormap('autumn');
    a(1, :) = [.8 .8 .8];
    if numel(unique(c)) == 1
        for kk = 1:size(a, 1)
            a(kk, :) = [.8 .8 .8];
        end
    end
    colormap(a);
end