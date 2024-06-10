%% Read data
raw_data = h5read('grid_01.hdf5','/data');


%% Clean up data
th_lo = 500;
th_hi = 4000;

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

th_roi_sum_lo = 120e3;
th_roi_sum_hi = 172e3;

grid_nstep = [49 49];
grid_step = 4000/55;
grid_step_offset = [0 0];

centroid_wdw_half_len = [35 35];

interactive_plot = false;

%% Assemble sets of data per chip, per hexa and per module
% Assemble chip sets
chip_data = pimega_chip_data(raw_data_filtered, px_array, chip_array);

% Assemble hexa sets and add gaps between chips with nominal spacing
hexa_data = pimega_hexa_data(chip_data, chip_array, chip_gap, hexa_array);

hexa_ref = cell(size(hexa_data));
hexa_x_ref = cell(size(hexa_data));
hexa_y_ref = cell(size(hexa_data));
for i=1:numel(hexa_data)
    hexa_ref{i} = repmat(i, size(hexa_data{i}));
    hexa_x_ref{i} = 1:size(hexa_data{i},2);
    hexa_y_ref{i} = 1:size(hexa_data{i},1);
end

% Assemble module sets and add gaps between hexas with nominal spacing
module_data = pimega_module_data(hexa_data, hexa_array, hexa_gap, module_array);
module_hexa_ref = pimega_module_data(hexa_ref, hexa_array, hexa_gap, module_array);
module_hexa_x_ref = pimega_module_data(hexa_x_ref, hexa_array, hexa_gap, module_array);
module_hexa_y_ref = pimega_module_data(hexa_y_ref, hexa_array, hexa_gap, module_array);

%% Assemble PIMEGA 540D detector by applying nominal rotations of modules and add gaps between modules with nominal spacing
detector_data = pimega_540d_data(module_data, module_gap);
detector_hexa_ref = pimega_540d_data(module_hexa_ref, module_gap);

%% Estimate beam centroids
x_centroid = nan(grid_nstep(1), grid_nstep(2));
y_centroid = nan(grid_nstep(1), grid_nstep(2));
centroid_hexa_ref = nan(grid_nstep(1), grid_nstep(2));

for i=1:grid_nstep(1)
    for j=1:grid_nstep(2)
        try
            x = round((-centroid_wdw_half_len(1):centroid_wdw_half_len(1)) + grid_step_offset(1) + (i-1)*grid_step);
            y = round((-centroid_wdw_half_len(1):centroid_wdw_half_len(2)) + grid_step_offset(2) + (j-1)*grid_step);
            
            roi = detector_data(y, x);
            roi_sum = sum(sum(roi, 'omitnan'));
            
            roi_hexa = detector_hexa_ref(y, x);
            roi_hexa(isnan(roi_hexa)) = [];
            
            % Check if data used to estimate centroids spans only one hexa
            if ~isempty(roi_hexa)
                roi_only_one_hexa = all(all(roi_hexa == roi_hexa(1)));
            else
                roi_only_one_hexa = false;
            end
            
            if interactive_plot
                figure(200);
                imshow(roi/th_hi)
                title(sprintf('px = %d, %d (step = %d, %d)', x(1), y(1),  i,j))
            end
            
            if roi_sum > th_roi_sum_lo && roi_sum < th_roi_sum_hi && roi_only_one_hexa
                x_centroid(i,j) = sum(x.*sum(roi, 1, 'omitnan'))/roi_sum;
                y_centroid(i,j) = sum(y'.*sum(roi, 2, 'omitnan'))/roi_sum;
                
                centroid_hexa_ref(i,j) = detector_hexa_ref(round(y_centroid(i,j)), round(x_centroid(i,j)));
            end
        catch
        end
    end
end

%% Fit
xy_meas = [x_centroid(:), y_centroid(:)];

% Artificially rotate data to compensante for experimental grid rotation
% grid_angle = 0.0045;
% grid_xoff = 0;
% grid_yoff = 0;
% Mrot = [cos(grid_angle) sin(grid_angle) 0; -sin(grid_angle) cos(grid_angle) 0; grid_xoff grid_yoff 1];
% xy_meas = Mrot\[xy_meas'; ones(1,numel(x_centroid))];
% xy_meas = xy_meas(1:2,:)';

[y_target, x_target] = meshgrid((0:grid_nstep(1)-1)*grid_step, (0:grid_nstep(2)-1)*grid_step);
xy_target = [x_target(:) y_target(:)];

[M, xy] = pimega_regression(xy_meas, xy_target, centroid_hexa_ref(:));

%% Map
% [map_x, map_y] = meshgrid(1:3578, 1:3578); 
% [x,y]=meshgrid(1:3578, 1:3578);
% 
% for i=1:3578
%     for j=1:3578
%         if all(~isnan(detector_hexa_ref(i,j)))
%             map = M{detector_hexa_ref(i,j)}*[j; i; 1];
%             map_x(i,j) = map(1);
%             map_y(i,j) = map(2);
%         end
%     end
% end
% %surf(map_x, map_y, detector_hexa_ref); view(0,90)



%% Check results
figure;
imshow(detector_data/th_hi); colormap(gca,'parula')
hold all
plot(x_centroid(:), y_centroid(:), 'or');

figure;
plot(xy.meas(:,1), xy.meas(:,2), '+r');
hold all
plot(xy.target(:,1), xy.target(:,2), 'ob');
clr = lines(length(xy.blocks));
%leg = cell(length(xy.blocks),1);
for i = 1:length(xy.blocks)
    plot(xy.corrected(xy.block_idx{i},1), xy.corrected(xy.block_idx{i},2), '.', 'Color', clr(i,:));
    %plot3(xy.corrected(xy.block_idx{i},1), xy.corrected(xy.block_idx{i},2), abs(xy.dist(xy.block_idx{i})), '.', 'Color', clr(i,:)); view(0,90);
    %leg{i} = sprintf('Corrected (HEXA # %d)', i);
end
title('Corrections by HEXA')
legend('Measured', 'Target', 'Corrected')
set(gca, 'YDir', 'reverse')
   
figure;
plot(abs(xy.dist(:)))
title('Error (distance between target and corrected positions)');