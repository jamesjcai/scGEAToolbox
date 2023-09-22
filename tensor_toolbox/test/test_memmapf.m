rng('default')
randData = rand([5, 5, 3], 'single');
fileID = fopen('mybinary1.bin', 'w');
fwrite(fileID, randData, 'single');
fclose(fileID);

%%

m = memmapfile('mybinary1.bin', ...
    'Format', {'single', [5, 5, 3], 'x'});
A = m.Data.x;

m.Writable = true;

X = uint16(1:1:15);
m.Data = X;

m.Offset = 0;
m.Repeat = 35;
%reshape(m.Data,5,7)';

% https://www.mathworks.com/help/matlab/import_export/writing-to-a-mapped-file.html