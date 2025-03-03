function [y] = e_runningasaddons

y1 = contains(which('scgeatool'),'MATLAB Add-Ons');
t = matlab.addons.installedAddons;
[y2,b]=ismember('scGEAToolbox', t.Name);
if y2
    y3 = t.Enabled(b);
else
    y3 = false;
end
y = y1 & y2 & y3;
