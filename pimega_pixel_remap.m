function [pxremap, hexa_centers, x_hexa, y_hexa] = pimega_pixel_remap(det)

nx = det.chip_array(2)*det.px_array(2);
ny = det.chip_array(1)*det.px_array(1);

% Mapping data
[ycoord, xcoord] = meshgrid(1:ny, 1:nx);
[module_data_x, hexa_data_x] = pimega_module_data_from_raw(xcoord, det);
pxremap.x = pimega_540d_data(module_data_x, det);                               % FIXME: should not be 540D-specifc
pxremap.y = pimega_540d_data(pimega_module_data_from_raw(ycoord, det), det);    % FIXME: should not be 540D-specifc

if nargout > 1
    % Reference number of each hexa
    hexa_ref = cell(hexa_data_x);
    for i=1:numel(hexa_data_x)
        hexa_ref{i} = repmat(i, size(hexa_data_x{i}));
    end
    detector_hexa_ref = pimega_540d_data(pimega_module_data(hexa_ref, det), det); % FIXME: should not be 540D-specifc
    [x_hexa,y_hexa,hexa_centers] = find_hexa_position(detector_hexa_ref, 1:numel(hexa_data_x));
end