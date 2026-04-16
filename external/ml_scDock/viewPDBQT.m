function viewPDBQT(filename)
    % viewPDBQT_withBonds - MATLAB PDBQT file viewer with bond inference
    % Usage: viewPDBQT_withBonds('molecule.pdbqt')

    %% --- Input Validation ---
    if nargin < 1 || ~ischar(filename)
        error('Please provide a valid PDBQT filename as a string.');
    end
    if ~isfile(filename)
        error('File "%s" not found.', filename);
    end

    %% --- Read PDBQT File ---
    fid = fopen(filename, 'r');
    if fid == -1
        error('Unable to open file: %s', filename);
    end

    atomCoords = [];
    atomTypes = {};

    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, 'ATOM') || startsWith(line, 'HETATM')
            % Parse coordinates (PDBQT format is fixed-width)
            x = str2double(line(31:38));
            y = str2double(line(39:46));
            z = str2double(line(47:54));
            atomType = strtrim(line(77:end)); % Atom type at end

            atomCoords(end+1, :) = [x, y, z];
            atomTypes{end+1} = atomType;
        end
    end
    fclose(fid);

    %% --- Assign Colors by Atom Type ---
    uniqueTypes = unique(atomTypes);
    cmap = lines(length(uniqueTypes));
    colors = zeros(length(atomTypes), 3);
    for i = 1:length(atomTypes)
        idx = find(strcmp(uniqueTypes, atomTypes{i}));
        colors(i, :) = cmap(idx, :);
    end

    %% --- Infer Bonds ---
    % Approximate covalent radii (Å) for common atoms
    covalentRadii = struct('H', 0.31, 'C', 0.76, 'N', 0.71, 'O', 0.66, ...
                           'P', 1.07, 'S', 1.05, 'Cl', 1.02, 'F', 0.57);
    bondPairs = [];
    nAtoms = size(atomCoords, 1);

    for i = 1:nAtoms-1
        for j = i+1:nAtoms
            type1 = atomTypes{i};
            type2 = atomTypes{j};

            % Extract element symbol (first letter(s))
            elem1 = regexp(type1, '^[A-Za-z]+', 'match', 'once');
            elem2 = regexp(type2, '^[A-Za-z]+', 'match', 'once');

            if isfield(covalentRadii, elem1) && isfield(covalentRadii, elem2)
                maxBondDist = covalentRadii.(elem1) + covalentRadii.(elem2) + 0.4; % tolerance
                dist = norm(atomCoords(i,:) - atomCoords(j,:));
                if dist <= maxBondDist
                    bondPairs(end+1, :) = [i, j];
                end
            end
        end
    end

    %% --- Plot Molecule ---
    figure('Name', sprintf('PDBQT Viewer with Bonds - %s', filename), 'Color', 'w');
    hold on;

    % Plot bonds first
    for k = 1:size(bondPairs, 1)
        i = bondPairs(k, 1);
        j = bondPairs(k, 2);
        plot3([atomCoords(i,1), atomCoords(j,1)], ...
              [atomCoords(i,2), atomCoords(j,2)], ...
              [atomCoords(i,3), atomCoords(j,3)], ...
              'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    end

    % Plot atoms
    scatter3(atomCoords(:,1), atomCoords(:,2), atomCoords(:,3), ...
             80, colors, 'filled', 'MarkerEdgeColor', 'k');

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('PDBQT Molecular Structure with Bonds');
    axis equal;
    grid on;
    rotate3d on;

    % Legend
    legend(uniqueTypes, 'Location', 'bestoutside');
    hold off;
end
