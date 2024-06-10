function [M, xy] = pimega_regression(xy_meas, xy_target, block_number)

if isempty(block_number) || nargin < 3
    block_number = ones(size(xy_meas,1));
end

A = [xy_meas xy_target block_number];

% Remove any row with NaN value
idx = isnan(A')';
idx = any(idx')';
A(idx,:) = [];

% Extend measured and target vectors with one to cope with translation
% (besides rotation) - Transpose to comply with typical transformation
% matrix convetion
npts = size(A,1);
xy_meas_ = [A(:,1:2)'; ones(1,npts)];
xy_target_ = [A(:,3:4)'; ones(1,npts)];

block_number_nonan = A(:,5);
blocks = unique(block_number_nonan);
nblocks = length(blocks);

xy_corrected_ = zeros(size(xy_meas_));
block_idx = cell(nblocks);
M = cell(nblocks,1);
for i=1:nblocks
    idx = block_number_nonan == blocks(i);
    M{i} = xy_meas_(:, idx)*pinv(xy_target_(:, idx));
    % M{i} = xy_target_(:, idx)*pinv(xy_meas_(:, idx)); % inverse matrix
    xy_corrected_(:, idx) = M{i}\xy_meas_(:, idx);
    block_idx{i} = idx;
end

xy.corrected = xy_corrected_(1:2,:)';
xy.meas = xy_meas_(1:2,:)';
xy.target = xy_target_(1:2,:)';
xy.dist = sqrt((xy.corrected(:,1)-xy.target(:,1)).^2 + (xy.corrected(:,2)-xy.target(:,2)).^2);
xy.blocks = blocks;
xy.block_idx = block_idx;