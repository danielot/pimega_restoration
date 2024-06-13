function chip_data = pimega_chip_data(raw_data, detector)

chip_data = cell(detector.chip_array);
for i=1:detector.chip_array(1)
    for j=1:detector.chip_array(2)
        chip_data{i,j} = raw_data((i-1)*detector.px_array(1) + (1:detector.px_array(1)), (j-1)*detector.px_array(2) + (1:detector.px_array(2)));
    end
end