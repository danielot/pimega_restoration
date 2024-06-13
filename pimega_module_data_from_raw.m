function module_data = pimega_module_data_from_raw(raw_data, detector, interp_method)

if nargin < 3
    interp_method = 'nearest';
end

chip_data = pimega_chip_data(raw_data, detector);
hexa_data = pimega_hexa_data(chip_data, detector);
if detector.hexa_tilt ~= 0
    hexa_data = regrid_hexa_data(hexa_data, detector, interp_method);
end
module_data = pimega_module_data(hexa_data, detector);


function chip_data = pimega_chip_data(raw_data, detector)

chip_data = cell(detector.chip_array);
for i=1:detector.chip_array(1)
    for j=1:detector.chip_array(2)
        chip_data{i,j} = raw_data((i-1)*detector.px_array(1) + (1:detector.px_array(1)), (j-1)*detector.px_array(2) + (1:detector.px_array(2)));
    end
end

function hexa_data = regrid_hexa_data(hexa_data, detector, interp_method)

xregrid = -(detector.px_array(1)-1)/2:(detector.px_array(1)-1)/2;
yregrid = -(5*detector.chip_gap + 6*detector.px_array(2)-1)/2:(5*detector.chip_gap + 6*detector.px_array(2)-1)/2;  % FIXME: harcoded values 5 an 6
[xrealmesh, yrealmesh] = meshgrid(xregrid*cos(detector.hexa_tilt),yregrid*detector);
[xregridmesh, yregridmesh] = meshgrid(xregrid,yregrid);
for i=1:numel(hexa_data)
    hexa_data{i} = griddata(xrealmesh, yrealmesh, double(hexa_data{i}), xregridmesh, yregridmesh, interp_method);
    hexa_data{i}(hexa_data{i}<0) = 0;
end

function module_data = pimega_module_data(hexa_data, detector)

nx_module = detector.hexa_array(1)/detector.module_array(1);
ny_module = detector.hexa_array(2)/detector.module_array(2);

if ~iscell(detector.hexa_gap) && isscalar(detector.hexa_gap)
    detector.hexa_gap_ = cell(detector.module_array);
    detector.hexa_gap_(:) = {repmat(detector.hexa_gap,1,ny_module-1)};
    detector.hexa_gap = detector.hexa_gap_;
end

module_data = cell(detector.module_array);
for i=1:detector.module_array(1)
    for j=1:detector.module_array(2)
        module = cell(1);
        module(1:2:2*ny_module-1) = hexa_data((i-1)*nx_module + (1:nx_module), (j-1)*ny_module + (1:ny_module));
        for k=1:ny_module-1
            module{2*k} = nan(size(hexa_data{1},1), detector.hexa_gap{i,j}(k));
        end
        module_data{i,j} = horzcat(module{:});
    end
end