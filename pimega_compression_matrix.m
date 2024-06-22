function M = pimega_compression_matrix(n, f, interp_method)
% M = PIMEGA_COMPRESSION_MATRIX(n,f)
%
%
% For an even n, M is:
%  n/2-(n/2-1)/f  ...  0         0      0       0         0        0    . . .   0
%        .       .     .         .      .       .         .        .            .
%        .        .    .         .      .       .         .        .            .
%        .         .   .         .      .       .         .        .            .
%        0           1/f-1       1      0       0         0        0            0
%        0             0         0      1     1/f-1       0        0            0
%        0    . . .    0         0      0     2-1/f     2/f-2      0    . . .   0
%        0             0         0      0       0       3-2/f    3/f-3          0
%        0             0         0      0       0         0      4-3/f          0
%        .             .         .      .       .         .        .   .        .
%        .             .         .      .       .         .        .    .       .
%        .             .         .      .       .         .        .     .      .
%        0    . . .    0         0      0       0         0        0  ...  n/2-(n/2-1)/f
%
% Sum across rows of M is 1.
% Sum across columns of M is 1/f.


if rem(n,2) == 0
    if n*f < n-2
        error('Only valid for small compressions. Maximum lost of two rows only (upper and lower most).');
    end
    invf = 1/f;
    seq = (0:n/2-1);
    seqinvf = seq*invf;
    M = diag(seq+1-seqinvf);
    M_ = zeros(n/2);
    M_(1:end-1, 2:end) = diag(seqinvf(2:end)-seq(2:end));
    M = M+M_;
    M = blkdiag(rot90(rot90(M)), M);
end

if nargin > 2 && strcmpi(interp_method, 'nearest')
    M = round(M*f);
end

M(all(M == 0, 2), :) = nan;