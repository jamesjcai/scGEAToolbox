function [totalRAM, availableRAM, usedRAM] = check_memory()
    % CHECK_MEMORY: Returns total, available, and used physical memory (RAM) across platforms.
    % Outputs:
    %   totalRAM     -> Total physical RAM (in GB)
    %   availableRAM -> Available physical RAM (in GB)
    %   usedRAM      -> Used physical RAM (in GB)
    % AUTHOR: Selim Romero, Texas A&M University

    if ispc
        % Windows: Use memory function
        [~, systemview] = memory;

        % Total system RAM (physical memory)
        totalRAM = systemview.PhysicalMemory.Total / 1e9; % Convert to GB

        % Available physical RAM
        availableRAM = systemview.PhysicalMemory.Available / 1e9; % Convert to GB

        % Used physical RAM
        usedRAM = totalRAM - availableRAM;

    elseif isunix
        if ismac
            % macOS: Use vm_stat
            [totalRAM, availableRAM, usedRAM] = mac_memory();
        else
            % Linux: Use free command
            [totalRAM, availableRAM, usedRAM] = linux_memory();
        end

    else
        error('Unsupported operating system.');
    end
end

function [totalRAM, availableRAM, usedRAM] = linux_memory()
    % Helper function to fetch physical memory stats on Linux
    [~, cmdout] = system('free -b'); % Get memory stats in bytes
    lines = strsplit(cmdout, '\n');
    memLine = strsplit(lines{2}); % Second line has memory stats

    % Parse values (in bytes)
    totalRAM = str2double(memLine{2}) / 1e9;     % Total physical memory
    availableRAM = str2double(memLine{7}) / 1e9; % Available physical memory
    usedRAM = totalRAM - availableRAM;          % Used physical memory
end

function [totalRAM, availableRAM, usedRAM] = mac_memory()
    % Helper function to fetch physical memory stats on macOS
    [~, cmdout] = system('vm_stat');
    lines = splitlines(cmdout);
    pageSize = 4096; % Default page size in bytes for macOS

    % Initialize values
    freePages = 0;
    activePages = 0;
    inactivePages = 0;
    wiredPages = 0;

    % Parse vm_stat output
    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if startsWith(line, 'Pages free:')
            freePages = sscanf(line, 'Pages free: %d');
        elseif startsWith(line, 'Pages active:')
            activePages = sscanf(line, 'Pages active: %d');
        elseif startsWith(line, 'Pages inactive:')
            inactivePages = sscanf(line, 'Pages inactive: %d');
        elseif startsWith(line, 'Pages wired down:')
            wiredPages = sscanf(line, 'Pages wired down: %d');
        end
    end

    % Calculate memory stats (in bytes)
    totalRAM = (freePages + activePages + inactivePages + wiredPages) * pageSize / 1e9; % Total physical memory
    usedRAM = (activePages + inactivePages + wiredPages) * pageSize / 1e9;             % Used physical memory
    availableRAM = freePages * pageSize / 1e9;                                        % Available physical memory
end
