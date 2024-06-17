dataset_filename = {
    'CAT/grid_module1_000.hdf5'
    'CAT/grid_module2_000.hdf5'
    'CAT/grid_module3_000.hdf5'
    'CAT/grid_module4_000.hdf5'
    %'CAT/grid_module1_center_000.hdf5'
    };

pixel_map_filename = {
    'CAT/pixel_remap_gap0.mat'
    %'CAT/pixel_remap_gap0_module3_2pxup.mat'
    %'CAT/pixel_remap_gap0_module3_1pxup_module4_1pxdown.mat'
    'CAT/pixel_remap_gap0_module1_1pxup_module3_1pxup_module4_1pxdown.mat'
    };


ndata = numel(dataset_filename);
npxremap = numel(pixel_map_filename);

pxremap_ = cell(npxremap, 1);
hexa_centers_ = cell(npxremap, 1);
x_hexa_ = cell(npxremap, 1);
y_hexa_ = cell(npxremap, 1);
for i=1:npxremap
    W = load(pixel_map_filename{i});
    pxremap_{i} = W.pxremap;
    hexa_centers_{i} = W.hexa_centers;
    x_hexa_{i} = W.x_hexa;
    y_hexa_{i} = W.y_hexa;
end

nplot = ndata*npxremap;
nploty = 3;
if nplot > nploty
    nplotx = ceil(nplot/nploty);
else
    nplotx = 1;
    nploty = nplot;
end
ax = zeros(1, nplot);
ax2 = zeros(1, nplot);
k = 1;

for data_i = 1:ndata
    %% Read data
    raw_data = h5read(sprintf('%s', dataset_filename{data_i}),'/entry/data/data');
    
    %% Detector Parameters
    det = detparam('pimega_540D', '2', '');
    
    %% Sample parameters
    sample_detector_distance = 7;           % [meters]
    dft_lobes_step = [13.14 13.11];       % [pixels]
    
    switch dataset_filename{data_i}
        case 'CAT/grid_module1_000.hdf5'
            grid_offset = [111.4 69.6];         % [pixels]
            grid_offset_hexa_ref = 7;
            grid_rot_angle = -0.00355;          % [radians]
            grid_shear_angle = -0.006;          % [radians]
            
        case 'CAT/grid_module2_000.hdf5'
            grid_offset = [-40.5 83.4];         % [pixels]
            grid_offset_hexa_ref = 6;
            grid_rot_angle = -0.00355;          % [radians]
            grid_shear_angle = -0.006;          % [radians]
            
        case 'CAT/grid_module3_000.hdf5'
            grid_offset = [-26.7 47.3];         % [pixels]
            grid_offset_hexa_ref = 19;
            grid_rot_angle = -0.00355;          % [radians]
            grid_shear_angle = -0.006;          % [radians]
            
        case 'CAT/grid_module4_000.hdf5'
            grid_offset = [78.2 -68.55];        % [pixels]
            grid_offset_hexa_ref = 18;
            grid_rot_angle = -0.00355;          % [radians]
            grid_shear_angle = -0.006;          % [radians]
            
        case 'CAT/grid_module1_center_000.hdf5'
            grid_offset = [617.9 -2.27];         % [pixels]
            grid_offset_hexa_ref = 11;
            grid_rot_angle = -0.00355;          % [radians]
            grid_shear_angle = -0.006;          % [radians]
    end
    
    for j = 1:npxremap
        
        %% Assemble sets of data per chip, per hexa and per module
        % Detector data
        pxremap = pxremap_{j};
        hexa_centers = hexa_centers_{j};
        x_hexa = x_hexa_{j};
        y_hexa = y_hexa_{j};
        
        detector_data = apply_restoration(raw_data, pxremap);
        
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
        xyrottrans_byhexa = cell(numel(hexa_centers), 1);
        for i=1:size(hexa_centers,1)
            filt_idx_x = xyrottrans_all(1,:) >= (1-f)*x_hexa(i,1) & xyrottrans_all(1,:) <= (1+f)*x_hexa(i,2);
            filt_idx_y = xyrottrans_all(2,:) >= (1-f)*y_hexa(i,1) & xyrottrans_all(2,:) <= (1+f)*y_hexa(i,2);
            xyrottrans_byhexa{i} = xyrottrans_all(:, filt_idx_x & filt_idx_y);
        end
        
        %% Check results
        
        detector_data(detector_data == -1) = nan;
        
        [ny, nx] = size(detector_data);
        
        nx_off = hexa_centers(grid_offset_hexa_ref,1) + grid_offset(1);
        ny_off = hexa_centers(grid_offset_hexa_ref,2) + grid_offset(2);
        
detector_data(isnan(detector_data)) = 0;
detector_data = imgaussfilt((double(detector_data)+1),2)+1;
        
        figure(2000);
        ax(k) = subplot(nplotx, nploty, k);
        imagesc((1:nx)-nx_off, (1:ny)-ny_off, log(double(detector_data)));
        colormap(gca, 'bone')
        hold all
        title(dataset_filename{data_i}, 'Interpreter', 'none')
        
        
        for i=1:size(hexa_centers,1)
            plot(xyrottrans_byhexa{i}(1,:)-nx_off, xyrottrans_byhexa{i}(2,:)-ny_off, '+', 'LineWidth', 2);
        end
        axis equal
        
        
        figure(2001);
        ax2(k) = subplot(nplotx, nploty, k);
        imagesc(log(double(detector_data)));
        colormap(gca, 'bone')
        hold all
        title(dataset_filename{data_i}, 'Interpreter', 'none')
        
        
        for i=1:size(hexa_centers,1)
            plot(xyrottrans_byhexa{i}(1,:), xyrottrans_byhexa{i}(2,:), '+', 'LineWidth', 2);
        end
        axis equal

        k = k+1;

        
        
    end
end

linkaxes(ax)
linkaxes(ax2)

n = 50;
axis(ax(1), [-n n -n n])
axis(ax2(1), [-n n -n n])


