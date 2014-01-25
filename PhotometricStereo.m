function [rho, n] = PhotometricStereo(I, mask, L)

% [rho, n] = PhotometricStereo(I, mask, L)
%
% Run photometric stereo with given images 'I', shadow mask 'mask', and
% calibrated lighting 'L'.
%
% INPUT:
%   I: N1xN2xM array, with each level I(:,:,i) the i-th intensity image.
%   mask: N1xN2xM boolean array, with each level mask(:,:,i) the shadow mask
%         of image I(:,:,i), 0 for pixel being in shadow and 1 otherwise.
%   L: 3xM array for calibrated lighting.
%
% OUTPUT:
%   rho: N1xN2 matrix for the estimated albedo map.
%   n: N1xN2x3 array for the estimated normal map.
%
%   Author: Ying Xiong.
%   Created: Jan 25, 2014.

[N1, N2, M] = size(I);

b = nan(N1, N2, 3);
for j1 = 1:N1
  ShowProgress(j1, N1, 10, 'PhotometricStereo', 1);
  for j2 = 1:N2
    Ij = squeeze(I(j1,j2,:));
    tag = squeeze(mask(j1,j2,:));
    b(j1,j2,:) = (L(:,tag)') \ Ij(tag);
  end
end

rho = sqrt(sum(b.^2, 3));
n = b ./ repmat(rho, [1 1 3]);
