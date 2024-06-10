%% Read data
dataset_name = 'grid_module1_center_000';
raw_data = h5read(sprintf('CAT/%s.hdf5', dataset_name),'/entry/data/data');


%% Clean up the data
th_lo = 0;
th_hi = Inf;

raw_data_filtered = raw_data;
raw_data_filtered(raw_data_filtered > th_hi) = 0;
raw_data_filtered(raw_data_filtered < th_lo) = 0;

%% Parameters
px_array = [256 256];                   % [pixels]
chip_array = [12 12];                   % [pixels]
hexa_array = [2 12];                    % [pixels]
module_array = [2 2];                   % [pixels]

chip_gap = 3;                           % [pixels]
module_gap_x = [                        % [pixels]
    0 8; ...                
    4 3];
module_gap_y = [                        % [pixels]
    2 0; ...                
    6 7];
hexa_gap = {};                          % [pixels]
hexa_gap{1,1} = [1 0 1 0 0];
hexa_gap{2,1} = [0 1 0 0 0];
hexa_gap{1,2} = [0 0 0 0 0];
hexa_gap{2,2} = [0 0 0 0 2];
hexa_overlap_angle = 6.75*pi/180;       % [radians]

th_roi_sum_lo = 120e3;                  % [counts]
th_roi_sum_hi = 172e3;                  % [counts]

sample_detector_distance = 7;           % [meters]
dft_lobes_step = [13.140 13.122];       % [pixels]
pixel_size = 55e-6;                     % [meters]


switch dataset_name
    case 'grid_module1_000'
        grid_range_x = [-67 168];           % [pixels]
        grid_range_y = [-73 163];           % [pixels]
        grid_offset = [111.4 69.6];         % [pixels]
        grid_offset_hexa_ref = 7;
        grid_rot_angle = -0.00355;          % [radians]
        grid_shear_angle = -0.006;          % [radians]
        
    case 'grid_module2_000'
        grid_range_x = [-200 80];           % [pixels]
        grid_range_y = [-70 181];           % [pixels]
        grid_offset = [-40.5 83.4];           % [pixels]
        grid_offset_hexa_ref = 6;
        grid_rot_angle = -0.00355;          % [radians]
        grid_shear_angle = -0.006;          % [radians]
        
    case 'grid_module3_000'
        grid_range_x = [-80 200];           % [pixels]
        grid_range_y = [-181 55];           % [pixels]
        grid_offset = [-26.7 47.3];         % [pixels]
        grid_offset_hexa_ref = 19;
        grid_rot_angle = -0.00355;          % [radians]
        grid_shear_angle = -0.006;          % [radians]
        
    case 'grid_module4_000'
        grid_range_x = [-200 60];           % [pixels]
        grid_range_y = [-200 70];           % [pixels]
        grid_offset = [78.2 -68.55];        % [pixels]
        grid_offset_hexa_ref = 18;
        grid_rot_angle = -0.00355;          % [radians]
        grid_shear_angle = -0.006;          % [radians]
        
    case 'grid_module1_center_000'
        grid_range_x = [-170 150];           % [pixels]
        grid_range_y = [-190 150];           % [pixels]
        grid_offset = [617.9 -2.27];   % [pixels]
        grid_offset_hexa_ref = 11;
        grid_rot_angle = -0.00355;          % [radians]
        grid_shear_angle = -0.006;          % [radians]
end

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
detector_data = pimega_540d_data(module_data, module_gap_x, module_gap_y);
detector_hexa_ref = pimega_540d_data(module_hexa_ref, module_gap_x, module_gap_y);
[x_hexa,y_hexa,hexa_centers] = find_hexa_position(detector_hexa_ref, 1:numel(hexa_data));

%%
% Take into account conical beam
angle_step = atan(dft_lobes_step*pixel_size/sample_detector_distance);
x = tan((grid_range_x(1):grid_range_x(end))*angle_step(1))*sample_detector_distance/pixel_size;
y = tan((grid_range_y(1):grid_range_y(end))*angle_step(2))*sample_detector_distance/pixel_size;

[xmesh, ymesh] = meshgrid(x,y);
xy = [xmesh(:)'; ymesh(:)'];
Mrot = [cos(grid_rot_angle) -sin(grid_rot_angle); sin(grid_rot_angle) cos(grid_rot_angle)];
Mhshear = [1 tan(grid_shear_angle); 0 1];
xyrot = Mrot*Mhshear*xy;
xyrottrans = xyrot + repmat((grid_offset+hexa_centers(grid_offset_hexa_ref,:))',1, size(xyrot,2));

f= 0.001;
xyrottrans_ = {};
for i=1:numel(hexa_data)
    aa = xyrottrans(1,:) >= (1-f)*x_hexa(i,1) & xyrottrans(1,:) <= (1+f)*x_hexa(i,2);
    bb = xyrottrans(2,:) >= (1-f)*y_hexa(i,1) & xyrottrans(2,:) <= (1+f)*y_hexa(i,2);
    xyrottrans_{i} = xyrottrans(:, aa & bb);
end

alpha = 1/cos(hexa_overlap_angle);
ddddd = [
    1 alpha
    alpha 1
    1 alpha
    alpha 1
    1 alpha
    alpha 1
    1 alpha
    alpha 1
    1 alpha
    alpha 1
    1 alpha
    alpha 1
    alpha 1
    1 alpha
    alpha 1
    1 alpha
    alpha 1
    1 alpha
    alpha 1
    1 alpha
    alpha 1
    1 alpha
    alpha 1
    1 alpha
    ];

%% Check results
figure;
imagesc(log(double(detector_data)));
colormap(gca, 'bone')
hold all

for i=1:numel(hexa_data)
    M = [ddddd(i,1) 0; 0 ddddd(i,2)];
    xyrottrans__  = M*(xyrottrans_{i} - hexa_centers(i,:)') + hexa_centers(i,:)';
    
    plot(xyrottrans__(1,:), xyrottrans__(2,:), 'o');
    axis equal
end