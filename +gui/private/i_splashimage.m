function splashpng = i_splashimage()
%I_SPLASHIMAGE Path of today's splash picture, or "" if there is none.
%   splashpng = i_splashimage() picks one file from assets/Images/splash_folder,
%   the same one all day: the pick is seeded by the date. Shared by
%   gui.sc_splashscreen and gui.sc_simplesplash.

splashdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', ...
    'assets', 'Images', 'splash_folder');
a = dir(splashdir);
a = a(~[a.isdir]); % keep files only
splashpng = "";
if isempty(a), return; end

d = datetime('today');
seed = year(d)*10000 + month(d)*100 + day(d);
% Save and restore the caller's random stream. Seeding the picture-of-the-day
% pick is fine; leaving the session parked on that seed is not -- it then
% governs every later tsne, umap and clustering call in the session.
rngState = rng();
restoreRng = onCleanup(@() rng(rngState));
rng(seed);
idx = randi(numel(a));
splashpng = string(fullfile(a(idx).folder, a(idx).name));
end
