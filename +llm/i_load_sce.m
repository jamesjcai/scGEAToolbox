function sce = i_load_sce(sample_id, data_dir)
% LLM.I_LOAD_SCE  Locate and load cleandata.mat for a GEO sample.
%
%   sce = llm.i_load_sce(sample_id, data_dir)
%
%   Searches data_dir/<any_study_id>/<sample_id>/cleandata.mat first,
%   then falls back to data_dir/<sample_id>/cleandata.mat.
%   Errors if the file cannot be found under either layout.

hits = dir(fullfile(data_dir, '*', sample_id, 'cleandata.mat'));
if isempty(hits)
    flat = fullfile(data_dir, sample_id, 'cleandata.mat');
    if isfile(flat)
        mat_path = flat;
    else
        error('llm:i_load_sce:fileNotFound', ...
            'Cannot find cleandata.mat for sample "%s" under "%s".', ...
            sample_id, data_dir);
    end
else
    mat_path = fullfile(hits(1).folder, hits(1).name);
end
fprintf('Loading %s\n', mat_path);
s = load(mat_path, 'sce');
sce = s.sce;
end
