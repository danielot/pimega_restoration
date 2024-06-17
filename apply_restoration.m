function data_restored = apply_restoration(raw_data, pixel_remap)

% nx = size(raw_data,1);
% ny = size(raw_data,2);
% nframes = size(raw_data, 3);
% nraw = nx*ny;
% 
% data_restored = repmat(-1, [size(pixel_remap.x) nframes]);
% 
% idx_raw = (pixel_remap.x(:)-1) + (pixel_remap.y(:)-1)*ny + 1;
% idx_good = ~isnan(idx_raw(:)) & idx_raw(:) > -1;
% idx_good = repmat(idx_good, nframes, 1);
% for i=1:nframes
%     data_restored(idx_good) = raw_data(idx_raw(idx_good) + (i-1)*nraw);
% end

data_restored = repmat(-1, size(pixel_remap.x));
idx_raw = (pixel_remap.x(:)-1) + (pixel_remap.y(:)-1)*size(raw_data,2) + 1;
idx_good = idx_raw(:) > 0;
data_restored(idx_good) = raw_data(idx_raw(idx_good));