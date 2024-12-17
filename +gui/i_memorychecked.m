function [continue_to_try, prepare_input_only] = i_memorychecked(ram_needed)


if nargin < 1, ram_needed = 64; end

continue_to_try = false;
prepare_input_only = false;

[~, sys] = memory;
totalMemoryGB = sys.PhysicalMemory.Total / (1024^3);
s = sprintf('At least %d GB of memory is recommended to run this function. Detected total physical memory: %.2f GB.', ...
    ram_needed, totalMemoryGB);

if totalMemoryGB < ram_needed * 0.95
    answer = questdlg(s, ...
        '', 'Continue, but Only Prepare Input', 'Proceed Anyway', 'Cancel', ...
        'Continue, but Only Prepare Input');
    switch answer
        case 'Continue, but Only Prepare Input'
            continue_to_try = true;
            prepare_input_only = true;
        case 'Proceed Anyway'
            continue_to_try = true;            
        case 'Cancel'
            return;
        otherwise
            return;
    end
else
    continue_to_try = true;
end
