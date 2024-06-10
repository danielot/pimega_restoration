function chip_data = pimega_chip_data(raw_data, px_array, chip_array)

chip_data = cell(chip_array);
for i=1:chip_array(1)
    for j=1:chip_array(2)
        chip_data{i,j} = raw_data((i-1)*px_array(1) + (1:px_array(1)), (j-1)*px_array(2) + (1:px_array(2)));
    end
end