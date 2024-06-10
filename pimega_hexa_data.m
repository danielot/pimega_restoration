function hexa_data = pimega_hexa_data(chip_data, chip_array, chip_gap, hexa_array)

nx_hexa = chip_array(1)/hexa_array(1);
ny_hexa = chip_array(2)/hexa_array(2);

if ~iscell(chip_gap) && isscalar(chip_gap)
    chip_gap_ = cell(hexa_array);
    chip_gap_(:) = {repmat(chip_gap,1,nx_hexa-1)};
    chip_gap = chip_gap_;
end

hexa_data = cell(hexa_array);
for i=1:hexa_array(1)
    for j=1:hexa_array(2)
        hexa = cell(1);
        hexa(1:2:2*nx_hexa-1) = chip_data((i-1)*nx_hexa + (1:nx_hexa), (j-1)*ny_hexa + (1:ny_hexa));
        for k=1:nx_hexa-1
            hexa{2*k} = nan(chip_gap{i,j}(k), size(chip_data{1},2));
        end
        hexa_data{i,j} = vertcat(hexa{:});
    end
end