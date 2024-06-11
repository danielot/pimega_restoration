function hexa_data = regrid_hexa_data(hexa_data, px_array, chip_gap, scaling)

xregrid = -(px_array(1)-1)/2:(px_array(1)-1)/2;
yregrid = -(5*chip_gap + 6*px_array(2)-1)/2:(5*chip_gap + 6*px_array(2)-1)/2;  % FIXME: harcoded values 5 an 6
[xrealmesh, yrealmesh] = meshgrid(xregrid*scaling(1),yregrid*scaling(2));
[xregridmesh, yregridmesh] = meshgrid(xregrid,yregrid);
for i=1:numel(hexa_data)
    hexa_data{i} = griddata(xrealmesh, yrealmesh, double(hexa_data{i}), xregridmesh, yregridmesh, 'nearest');
    hexa_data{i}(hexa_data{i}<0) = 0;
end