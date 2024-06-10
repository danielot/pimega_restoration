function [x,y,c] = find_hexa_position(detector_hexa_ref, hexa_labels)

% Only for hexas aligned with X and Y axes
n = length(hexa_labels);
x = zeros(n,2);
y = zeros(n,2);
for i=1:n
    A = detector_hexa_ref == hexa_labels(i);
    x_ = find(any(A));
    y_ = find(any(A'));
    x(i,:) = [min(x_) max(x_)];
    y(i,:) = [min(y_) max(y_)];
end

c = [x(:,1)+x(:,2) y(:,1)+y(:,2)]/2;