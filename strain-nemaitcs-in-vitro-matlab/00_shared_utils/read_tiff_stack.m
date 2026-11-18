function stack = read_tiff_stack(fname, frames)
% READ_TIFF_STACK  Read a single- or multi-page TIFF into an H x W x N array.
%   stack = read_tiff_stack(fname)          all frames
%   stack = read_tiff_stack(fname, 1:100)   selected frames
% Keeps the native class (uint8 / uint16 / single ...). Uses imfinfo/imread
% only, so it works for any bit depth in MATLAB and GNU Octave.
% (loadtiff.m in this folder gives the same pixel values; it is kept because
% the historical scripts call it.)
info = imfinfo(fname);
n = numel(info);
if nargin < 2 || isempty(frames), frames = 1:n; end
if any(frames > n) || any(frames < 1)
    error('read_tiff_stack:range', '%s has %d frames; requested %d-%d.', ...
          fname, n, min(frames), max(frames));
end
im1 = imread(fname, frames(1));
if size(im1, 3) > 1, im1 = im1(:,:,1); end   % RGB -> first channel
stack = zeros(size(im1,1), size(im1,2), numel(frames), class(im1));
stack(:,:,1) = im1;
for k = 2:numel(frames)
    im = imread(fname, frames(k));
    if size(im, 3) > 1, im = im(:,:,1); end
    stack(:,:,k) = im;
end
end
