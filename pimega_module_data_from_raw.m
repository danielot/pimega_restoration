function [module_data, hexa_data, chip_data] = pimega_module_data_from_raw(raw_data, detector, interp_method)

if nargin < 3
    interp_method = 'nearest';
end

chip_data = pimega_chip_data(raw_data, detector);
hexa_data = pimega_hexa_data(chip_data, detector);
if detector.hexa_tilt ~= 0
    hexa_data = pimega_regrid_hexa_data(hexa_data, detector, interp_method);
end
module_data = pimega_module_data(hexa_data, detector);


