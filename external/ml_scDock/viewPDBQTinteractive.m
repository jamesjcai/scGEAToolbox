function viewPDBQTinteractive(filename)
    % interactivePDBQTViewer - Fully interactive PDBQT molecular viewer
    % Usage: interactivePDBQTViewer('molecule.pdbqt')

    %% --- Input Validation ---
    if nargin < 1 || ~ischar(filename)
        error('Please provide a valid PDBQT filename as a string.');
    end
    if ~isfile(filename)
        error('File "%s" not found.', filename);
    end

    %% --- Read PDBQT File ---
    [atomCoords, atomTypes] = readPDBQT(filename);

    %% --- Assign Colors by Atom Type ---
    uniqueTypes = unique(atomTypes);
    cmap = lines(length(uniqueTypes));
    colors = zeros(length(atomTypes), 3);
    for i = 1:length(atomTypes)
        idx = find(strcmp(uniqueTypes, atomTypes{i}));
        colors(i, :) = cmap(idx, :);
    end

    %% --- Infer Bonds ---
    bondPairs = inferBonds(atomCoords, atomTypes);

    %% --- Plot Molecule ---
    fig = figure('Name', sprintf('Interactive PDBQT Viewer - %s', filename), ...
                 'Color', 'w', 'NumberTitle', 'off');
    hold on;

    % Plot bonds
    for k = 1:size(bondPairs, 1)
        i = bondPairs(k, 1);
        j = bondPairs(k, 2);
        plot3([atomCoords(i,1), atomCoords(j,1)], ...
              [atomCoords(i,2), atomCoords(j,2)], ...
              [atomCoords(i,3), atomCoords(j,3)], ...
              'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    end

    % Plot atoms
    scatterHandle = scatter3(atomCoords(:,1), atomCoords(:,2), atomCoords(:,3), ...
                              80, colors, 'filled', 'MarkerEdgeColor', 'k');

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Interactive PDBQT Molecular Structure');
    axis equal;
    grid on;
    rotate3d on;
    legend(uniqueTypes, 'Location', 'bestoutside');

    %% --- Interaction State ---
    selectedAtoms = []; % store clicked atoms for measurement

    %% --- Click Callback ---
    set(fig, 'WindowButtonDownFcn', @onClick);

    function onClick(~, ~)
        % Get click location in 3D
        clickPoint = get(gca, 'CurrentPoint');
        clickPoint = clickPoint(1, 1:3);

        % Find nearest atom
        distances = vecnorm(atomCoords - clickPoint, 2, 2);
        [minDist, idx] = min(distances);

        if minDist < 1.0 % threshold for selection
            selectedAtoms(end+1) = idx;

            % Show atom info
            fprintf('Atom %d: Type=%s, Coordinates=(%.3f, %.3f, %.3f)\n', ...
                idx, atomTypes{idx}, atomCoords(idx,1), atomCoords(idx,2), atomCoords(idx,3));

            % If two atoms selected, measure distance
            if numel(selectedAtoms) == 2
                a1 = selectedAtoms(1);
                a2 = selectedAtoms(2);
                dist = norm(atomCoords(a1,:) - atomCoords(a2,:));
                fprintf('Distance between Atom %d and Atom %d: %.3f Å\n', a1, a2, dist);

                % Draw temporary measurement line
                plot3([atomCoords(a1,1), atomCoords(a2,1)], ...
                      [atomCoords(a1,2), atomCoords(a2,2)], ...
                      [atomCoords(a1,3), atomCoords(a2,3)], ...
                      'r--', 'LineWidth', 2);

                % Reset selection
                selectedAtoms = [];
            end
        end
    end
end

%% --- Helper: Read PDBQT ---
function [coords, types] = readPDBQT(filename)
    fid = fopen(filename, 'r');
    coords = [];
    types = {};
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, 'ATOM') || startsWith(line, 'HETATM')
            x = str2double(line(31:38));
            y = str2double(line(39:46));
            z = str2double(line(47:54));
            atomType = strtrim(line(77:end));
            coords(end+1, :) = [x, y, z];
            types{end+1} = atomType;
        end
    end
    fclose(fid);
end

%% --- Helper: Infer Bonds ---
function bondPairs = inferBonds(coords, types)
    covalentRadii = struct('H', 0.31, 'C', 0.76, 'N', 0.71, 'O', 0.66, ...
                           'P', 1.07, 'S', 1.05, 'Cl', 1.02, 'F', 0.57);
    bondPairs = [];
    nAtoms = size(coords, 1);
    for i = 1:nAtoms-1
        for j = i+1:nAtoms
            elem1 = regexp(types{i}, '^[A-Za-z]+', 'match', 'once');
            elem2 = regexp(types{j}, '^[A-Za-z]+', 'match', 'once');
            if isfield(covalentRadii, elem1) && isfield(covalentRadii, elem2)
                maxBondDist = covalentRadii.(elem1) + covalentRadii.(elem2) + 0.4;
                dist = norm(coords(i,:) - coords(j,:));
                if dist <= maxBondDist
                    bondPairs(end+1, :) = [i, j];
                end
            end
        end
    end
end
