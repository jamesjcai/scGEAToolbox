function i_setautumncolor(c,cmapname)

if nargin<2, cmapname='autumn'; end
    a = colormap(cmapname);
    a(1, :) = [.8 .8 .8];
    if numel(unique(c)) == 1
        for kk = 1:size(a, 1)
            a(kk, :) = [.8 .8 .8];
        end
    end
    colormap(a);
end