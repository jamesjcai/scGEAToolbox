% % Example data to save
data = rand(100, 10); % Example data, replace with your own
% 
% % Define file name
% filename = 'example_data.h5';
% 
% % Write data to HDF5 file
% h5create(filename, '/data', size(data));
% h5write(filename, '/data', data);
% 
% disp('Data saved successfully.');

% Define filename and dataset name
your_data = data;

filename = 'my_data.h5ad';
dataset_name = '/my_dataset';

% Create the file and dataset with appropriate data type
file_id = H5F.create(filename, 'H5F_CREATE', 'H5P_DEFAULT');
data_type_id = H5T.copy('H5T_NATIVE_DOUBLE'); % Change data type as needed
dataspace_id = H5S.create_simple(size(your_data)); % Replace with data dimensions
dataset_id = H5D.create(file_id, dataset_name, data_type_id, dataspace_id, 'H5P_DEFAULT');

% Close the data space and data type identifiers
H5S.close(dataspace_id);
H5T.close(data_type_id);
% Write your data to the dataset
H5D.write(dataset_id, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', your_data);

% Close the dataset and file identifiers
H5D.close(dataset_id);
H5F.close(file_id);
