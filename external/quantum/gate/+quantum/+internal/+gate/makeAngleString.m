function angleString = makeAngleString(angle, piStr)

%   Copyright 2021-2022 The MathWorks, Inc.

if angle == 0
    angleString = "0";
    return
end

[num,denom] = rat(angle/pi);

if denom <= 64
    % simplify the angle
    if num==1
        num = "";
    elseif num==-1
        num = "-";
    end
    if denom==1
        denom = "";
    else
        denom = "/"+string(denom);
    end
    angleString = string(num)+piStr+string(denom);
else
    angleString = string(angle);
end
