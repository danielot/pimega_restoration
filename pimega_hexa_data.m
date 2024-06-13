function hexa_data = pimega_hexa_data(chip_data, detector)

nx_hexa = detector.chip_array(1)/detector.hexa_array(1);
ny_hexa = detector.chip_array(2)/detector.hexa_array(2);

if ~iscell(detector.chip_gap) && isscalar(detector.chip_gap)
    chip_gap = cell(detector.hexa_array);
    chip_gap(:) = {repmat(detector.chip_gap,1,nx_hexa-1)};
    detector.chip_gap = chip_gap;
end

hexa_data = cell(detector.hexa_array);
for i=1:detector.hexa_array(1)
    for j=1:detector.hexa_array(2)
        hexa = cell(1);
        hexa(1:2:2*nx_hexa-1) = chip_data((i-1)*nx_hexa + (1:nx_hexa), (j-1)*ny_hexa + (1:ny_hexa));
        for k=1:nx_hexa-1
            hexa{2*k} = nan(chip_gap{i,j}(k), size(chip_data{1},2));
        end
        hexa_data{i,j} = vertcat(hexa{:});
    end
end