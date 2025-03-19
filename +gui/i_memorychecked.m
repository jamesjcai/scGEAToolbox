function [continue_to_try, prepare_input_only] = i_memorychecked(ram_needed, parentfig)

if nargin < 2, parentfig = []; end
if nargin < 1, ram_needed = 32; end

continue_to_try = false;
prepare_input_only = false;

    if ispc    
        [~, sys] = memory;
        totalMemoryGB = sys.PhysicalMemory.Total / (1024^3);
        s = sprintf(['%d GB of memory is recommended to run this function. ' ...
            'Detected total physical memory: %.2f GB.'], ...
            ram_needed, totalMemoryGB);
        
        if totalMemoryGB < ram_needed * 0.95
            answer = gui.myQuestdlg(parentfig, s, ...
                '', {'Continue, but Only Prepare Input', ...
                'Proceed Anyway', 'Cancel'}, ...
                'Continue, but Only Prepare Input');
            if isempty(answer), return; end
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
    else

        s = '32 GB of memory is recommended to run this function.';
        
        answer = gui.myQuestdlg(parentfig, s, ...
            '', {'Continue, but Only Prepare Input', ...
            'Proceed Anyway', 'Cancel'}, ...
            'Continue, but Only Prepare Input');
        if isempty(answer), return; end
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
    end

end