function hexa_data = regrid_hexa_data(hexa_data, detector, interp_method)

if nargin < 3
    interp_method = 'nearest';
end

nx_hexa = detector.chip_array(1)/detector.hexa_array(1);

xregrid = -(detector.px_array(1)-1)/2:(detector.px_array(1)-1)/2;
yregrid = -((nx_hexa-1)*detector.chip_gap + nx_hexa*detector.px_array(2)-1)/2:((nx_hexa-1)*detector.chip_gap + nx_hexa*detector.px_array(2)-1)/2;
[xrealmesh, yrealmesh] = meshgrid(xregrid*cos(detector.hexa_tilt),yregrid);
[xregridmesh, yregridmesh] = meshgrid(xregrid,yregrid);
for i=1:numel(hexa_data)
    hexa_data{i} = griddata(xrealmesh, yrealmesh, double(hexa_data{i}), xregridmesh, yregridmesh, interp_method);
    hexa_data{i}(hexa_data{i}<0) = 0;
end