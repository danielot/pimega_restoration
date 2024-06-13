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