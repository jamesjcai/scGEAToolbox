function [raw_data] = reverse_processing( dataname,  processed_data )
%REVERSE_PROCESSING Summary of this function goes here
%   Detailed explanation goes here

normed_data=2.^processed_data-1;

if (strcmp(dataname,'Mousebrain'))
        nlabels=7;
        load(['Temp/' dataname '_mcimputeProcessing_libsize_' num2str(nlabels) 'labels.mat']);
else
        load(['Temp/' dataname '_mcimputeProcessing_libsize.mat']);
end

raw_data=bsxfun(@times, normed_data, libsize) / median(libsize);
end

