function detector_data = pimega_540d_data(module_data, module_gap_x, module_gap_y)

if nargin < 3
    module_gap_y = module_gap_x;
end

if isscalar(module_gap_x)
    module_gap_x = [module_gap_x 0; module_gap_x 0];
end

if isscalar(module_gap_y)
    module_gap_y = [module_gap_y, 0; 0, module_gap_y];
end

mnmax = 0;
for i=1:numel(module_data)
    [n,m] = size(module_data{i});    
    if max(m,n) > mnmax
        mnmax = max(m,n);
    end
end

module1 = nan(mnmax,mnmax);
module2 = nan(mnmax,mnmax);
module3 = nan(mnmax,mnmax);
module4 = nan(mnmax,mnmax);

module1(mnmax-size(module_data{1,1},1)+1:end, mnmax-size(module_data{1,1},2)+1:end) = module_data{1,1};
module2(mnmax-size(module_data{2,1},1)+1:end, mnmax-size(module_data{2,1},2)+1:end) = module_data{2,1};
module3(mnmax-size(module_data{1,2},1)+1:end, mnmax-size(module_data{1,2},2)+1:end) = module_data{1,2};
module4(mnmax-size(module_data{2,2},1)+1:end, mnmax-size(module_data{2,2},2)+1:end) = module_data{2,2};

module1 = module1';
module2 = module2(:,end:-1:1);
module3 = module3(end:-1:1,:);
module4 = module4(end:-1:1,end:-1:1)';

mn = size(module1,1);

if module_gap_x(2,1) > module_gap_x(1,1)
    module1 = [nan(size(module1,1), module_gap_x(2,1)-module_gap_x(1,1)) module1 nan(size(module1,1), module_gap_x(1,1))];
    module3 = [module3 nan(size(module3,1), module_gap_x(2,1))];
else
    module1 = [module1 nan(size(module1,1), module_gap_x(1,1))];
    module3 = [nan(size(module3,1), module_gap_x(1,1)-module_gap_x(2,1)) module3 nan(size(module3,1), module_gap_x(2,1))];
end

if module_gap_x(2,2) > module_gap_x(1,2)
    module2 = [nan(size(module2,1), module_gap_x(1,2)) module2 nan(size(module2,1), module_gap_x(2,2)-module_gap_x(1,2))];
    module4 = [nan(size(module4,1), module_gap_x(2,2)) module4];
else
    module2 = [nan(size(module2,1), module_gap_x(1,2)) module2];
    module4 = [nan(size(module4,1), module_gap_x(2,2)) module4 nan(size(module4,1), module_gap_x(1,2)-module_gap_x(2,2))];
end

if module_gap_y(1,2) > module_gap_y(1,1)
    module1 = [nan(module_gap_y(1,2)-module_gap_y(1,1), size(module1,2)); module1; nan(module_gap_y(1,1), size(module1,2))];
    module2 = [module2; nan(module_gap_y(1,2), size(module2,2))];
else
    module1 = [module1; nan(module_gap_y(1,1), size(module1,2))];
    module2 = [nan(module_gap_y(1,1)-module_gap_y(1,2), size(module2,2)); module2; nan(module_gap_y(1,2), size(module2,2))];
end

if module_gap_y(2,2) > module_gap_y(2,1)
    module3 = [nan(module_gap_y(2,1), size(module3,2)); module3; nan(module_gap_y(2,2)-module_gap_y(2,1), size(module3,2))];
    module4 = [nan(module_gap_y(2,2), size(module4,2)); module4];
else
    module3 = [nan(module_gap_y(2,1), size(module3,2)); module3];
    module4 = [nan(module_gap_y(2,2), size(module4,2)); module4; nan(module_gap_y(2,1)-module_gap_y(2,2), size(module4,2))];
end

detector_data = [module1 module2; module3 module4];