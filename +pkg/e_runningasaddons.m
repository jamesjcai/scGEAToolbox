function [y] = e_runningasaddons

% Initialize output variable
y = false;
if ~(ismcc || isdeployed)
    y1 = contains(which('scgeatool'),'MATLAB Add-Ons');
    %#exclude matlab.addons.installedAddons
    t = matlab.addons.installedAddons;
    [y2,b]=ismember('scGEAToolbox', t.Name);
    if y2
        y3 = t.Enabled(b);
    else
        y3 = false;
    end
    y = y1 & y2 & y3;
end