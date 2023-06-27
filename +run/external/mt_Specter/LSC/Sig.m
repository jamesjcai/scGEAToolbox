% Author: Modified package by Van Hoan Do
function f = Sig(value, x, k)
    % ML based learn parameter 
    f = value + 1.0;
    if k <= 2
        f = f -1;
    end
    if 120 < x && x < 130 
        f = f-1;
    end 
    if 260 < x && x < 270 
        f = f-1;
    end
    if 1550 < x && x < 1605
        f = f -1;
    end 
    if 1870 < x && x < 1890 
        f = f - 1;
    end 
    if 2120 < x && x < 2135 
        f = f - 1;
    end 
    if 2610 < x && x < 1801
        f = f -1;
    end
    if 3001 < x && x < 3010 
        f = f - 1;
    end 
    if x > 14000 
        f = f -1;
    end
end
