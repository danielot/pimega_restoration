function hexa_data = pimega_hexa_compensate_tilt(hexa_data, detector, interp_method)
 
if nargin < 3
    interp_method = 'sum_neighbors';
end

M = pimega_compression_matrix(size(hexa_data{1}, 2), cos(detector.hexa_tilt), interp_method);

for i=1:numel(hexa_data)
     hexa_data{i} = hexa_data{i}*M';
end