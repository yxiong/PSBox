function PSFindLightDirection(topDir)

%   Author: Ying Xiong.
%   Created: Jan 24, 2014.

for iProbe = 1:2
  % Find all light probe images.
  probeDir = fullfile(topDir, ['LightProbe-' num2str(iProbe)]);
  imgFiles = dir(fullfile(probeDir, 'Image_*.JPG'));
  % Load circle data.
  circle = textread(fullfile(probeDir, 'circle_data.txt'));
  threshold = 250;
  % Process each image.
  nImgs = length(imgFiles);
  L = zeros(3, nImgs);
  figure;
  [nRows, nCols] = NumSubplotRowsColsFromTotal(nImgs);
  for iImg = 1:nImgs
    I = imread(fullfile(probeDir, imgFiles(iImg).name));
    I = I(end:-1:1, :, :);
    opts = struct('Visualize', 'on');
    subplot(nRows, nCols, iImg);
    L(:,iImg) = FindLightDirectionFromChromeSphere(I, circle, threshold, opts);
  end
  drawnow;
  % Write result to output.
  dlmwrite(fullfile(probeDir, 'light_directions.txt'), L, ...
           'delimiter', ' ', 'precision', '%10.6f');
end
