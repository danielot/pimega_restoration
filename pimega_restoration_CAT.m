%% Read data
dataset_name = 'grid_module1_center_000';
raw_data = h5read(sprintf('CAT/%s.hdf5', dataset_name),'/entry/data/data');

%% Detector Parameters
det = detparam('pimega_540D', '2', '');

%% Sample parameters
sample_detector_distance = 7;           % [meters]
dft_lobes_step = [13.140 13.122];       % [pixels]

switch dataset_name
    case 'grid_module1_000'
        grid_offset = [111.4 69.6];         % [pixels]
        grid_offset_hexa_ref = 7;
        grid_rot_angle = -0.00355;          % [radians]
        grid_shear_angle = -0.006;          % [radians]
        
    case 'grid_module2_000'
        grid_offset = [-40.5 83.4];         % [pixels]
        grid_offset_hexa_ref = 6;
        grid_rot_angle = -0.00355;          % [radians]
        grid_shear_angle = -0.006;          % [radians]
        
    case 'grid_module3_000'
        grid_offset = [-26.7 47.3];         % [pixels]
        grid_offset_hexa_ref = 19;
        grid_rot_angle = -0.00355;          % [radians]
        grid_shear_angle = -0.006;          % [radians]
        
    case 'grid_module4_000'
        grid_offset = [78.2 -68.55];        % [pixels]
        grid_offset_hexa_ref = 18;
        grid_rot_angle = -0.00355;          % [radians]
        grid_shear_angle = -0.006;          % [radians]
        
    case 'grid_module1_center_000'
        grid_offset = [617.9 -2.27];         % [pixels]
        grid_offset_hexa_ref = 11;
        grid_rot_angle = -0.00355;          % [radians]
        grid_shear_angle = -0.006;          % [radians]
end

%% Assemble sets of data per chip, per hexa and per module

% Mapping data 
[ycoord, xcoord] = meshgrid(1:size(raw_data, 2), 1:size(raw_data, 1));
[module_data_x, hexa_data_x] = pimega_module_data_from_raw(xcoord, det);
detector_data_x = pimega_540d_data(module_data_x, det);
module_data_y = pimega_module_data_from_raw(ycoord, det);
detector_data_y = pimega_540d_data(module_data_y, det);

% Detector data
module_data = pimega_module_data_from_raw(raw_data, det);
detector_data = pimega_540d_data(module_data, det);

% Reference number of each hexa
hexa_ref = cell(hexa_data_x);
for i=1:numel(hexa_data_x)
    hexa_ref{i} = repmat(i, size(hexa_data_x{i}));
end
module_hexa_ref = pimega_module_data(hexa_ref, det);
detector_hexa_ref = pimega_540d_data(module_hexa_ref, det);
[x_hexa,y_hexa,hexa_centers] = find_hexa_position(detector_hexa_ref, 1:numel(hexa_data_x));

%% Take into account conical beam
grid_range_x = ([1 size(detector_data,1)] - hexa_centers(grid_offset_hexa_ref,1) - grid_offset(1))/dft_lobes_step(1);
grid_range_x = [floor(grid_range_x(1)) ceil(grid_range_x(2))];
grid_range_y = ([1 size(detector_data,2)] - hexa_centers(grid_offset_hexa_ref,2) - grid_offset(2))/dft_lobes_step(2);
grid_range_y = [floor(grid_range_y(1)) ceil(grid_range_y(2))];

angle_step = atan(dft_lobes_step*det.pixel_size/sample_detector_distance);
x = tan((grid_range_x(1):grid_range_x(end))*angle_step(1))*sample_detector_distance/det.pixel_size;
y = tan((grid_range_y(1):grid_range_y(end))*angle_step(2))*sample_detector_distance/det.pixel_size;

[xmesh, ymesh] = meshgrid(x,y);
xy = [xmesh(:)'; ymesh(:)'];
Mrot = [cos(grid_rot_angle) -sin(grid_rot_angle); sin(grid_rot_angle) cos(grid_rot_angle)];
Mhshear = [1 tan(grid_shear_angle); 0 1];
xyrot = Mrot*Mhshear*xy;
xyrottrans_all = xyrot + repmat((grid_offset+hexa_centers(grid_offset_hexa_ref,:))',1, size(xyrot,2));

f= 0.001;
xyrottrans_byhexa = cell(numel(hexa_data_x), 1);
for i=1:numel(hexa_data_x)
    filt_idx_x = xyrottrans_all(1,:) >= (1-f)*x_hexa(i,1) & xyrottrans_all(1,:) <= (1+f)*x_hexa(i,2);
    filt_idx_y = xyrottrans_all(2,:) >= (1-f)*y_hexa(i,1) & xyrottrans_all(2,:) <= (1+f)*y_hexa(i,2);
    xyrottrans_byhexa{i} = xyrottrans_all(:, filt_idx_x & filt_idx_y);
end

%% Check results
figure;
imagesc(log(double(detector_data)));
colormap(gca, 'bone')
hold all

for i=1:numel(hexa_data_x)
    plot(xyrottrans_byhexa{i}(1,:), xyrottrans_byhexa{i}(2,:), '+', 'LineWidth', 2);    
end
axis equal
