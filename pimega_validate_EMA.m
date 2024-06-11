%% Read data
raw_data = double(h5read('EMA/raw_center.h5','/entry/data/data'));


%% Clean up data
th_lo = 0;
th_hi = 100000;

raw_data_filtered = raw_data;
raw_data_filtered(raw_data_filtered > th_hi) = 0;
raw_data_filtered(raw_data_filtered < th_lo) = 0;

%% Parameters
px_array = [256 256];
chip_array = [12 12];
hexa_array = [2 12];
module_array = [2 2];

chip_gap = 3;
hexa_gap = 50;
module_gap = 6;

%% Assemble sets of data per chip, per hexa and per module
chip_data = pimega_chip_data(raw_data_filtered, px_array, chip_array);
hexa_data = pimega_hexa_data(chip_data, chip_array, chip_gap, hexa_array);
module_data = pimega_module_data(hexa_data, hexa_array, hexa_gap, module_array);
detector_data = pimega_540d_data(module_data, module_gap);

%% Plot
[x,y] = meshgrid(1:3578, 1:3578);
%figure;
%imshow(griddata(map_x,map_y,detector_data_,x,y,'linear')/4000)
%colormap(gca, 'parula')

B=uint32(zeros(3578,3578));
for i=1:length(M)
    A = uint32(detector_data/4000*2^32-1);
    A(detector_hexa_ref ~= xy.blocks(i)) = 0;
    MM = (M{i});
    MM(3,:) = [0 0 1];
    MM = MM';
    xx = imwarp(A, affine2d(MM));
    dx=min(size(B,1), size(xx,1));
    dy=min(size(B,2), size(xx,2));
    
    B(1:dx,1:dy) = B(1:dx,1:dy) + xx(1:dx,1:dy);
    %imshow(A); colormap(gca, 'parula')
    %pause
    %
    %pause
end
imshow(B); colormap(gca, 'parula')