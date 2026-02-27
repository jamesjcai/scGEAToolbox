function obj = estimatepotency(obj, speciesid, forced)
if nargin < 3, forced = false; end
if nargin < 2, speciesid = []; end
if forced || ~obj.hasCellAttribute('cell_potency')
    if isempty(speciesid)
        FigureHandle=[];
        speciesid = gui.i_selectspecies(2, false, FigureHandle);

    elseif ischar(speciesid)
        [y, idx] = ismember(lower(speciesid), {'human', 'mouse'});
        if y, speciesid = idx; end
    elseif isstring(speciesid)
        [y, idx] = ismember(lower(speciesid), ["human", "mouse"]);
        if y, speciesid = idx; end
    elseif isnumeric(speciesid)
        if ~ismember(speciesid, [1, 2])
            error('SPECIESID should be 1 (human) or 2 (mouse)');
        end
    end
    r = sc_potency(obj.X, obj.g, speciesid);
    obj.setCellAttribute('cell_potency', r);
    disp('cell_potency added to SCE.LIST_CELL_ATTRIBUTES.');
else
    disp('cell_potency existed in SCE.LIST_CELL_ATTRIBUTES.');
end
end
