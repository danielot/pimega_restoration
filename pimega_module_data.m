function module_data = pimega_module_data(hexa_data, hexa_array, hexa_gap, module_array)

nx_module = hexa_array(1)/module_array(1);
ny_module = hexa_array(2)/module_array(2);

if ~iscell(hexa_gap) && isscalar(hexa_gap)
    hexa_gap_ = cell(module_array);
    hexa_gap_(:) = {repmat(hexa_gap,1,ny_module-1)};
    hexa_gap = hexa_gap_;
end

module_data = cell(module_array);
for i=1:module_array(1)
    for j=1:module_array(2)
        module = cell(1);
        module(1:2:2*ny_module-1) = hexa_data((i-1)*nx_module + (1:nx_module), (j-1)*ny_module + (1:ny_module));
        for k=1:ny_module-1
            module{2*k} = nan(size(hexa_data{1},1), hexa_gap{i,j}(k));
        end
        module_data{i,j} = horzcat(module{:});
    end
end