function snap = vec2snap(vec, xsize, ysize)
% vec2snap converts a reshaped vector snapshot into a two-d snap image
% 
% syntax: snap = vec2snap(vec, xsize, ysize)
% 
% input: vec - (xsize*ysize)-by-1 array, a spatial snapshot vector reshaped
%              from a two-d snapshot of scalar field.
%        xsize - int, size of the snap shot matrix in x direction
%        ysize - int, size of the snap shot matrix in y direction
%
% output: snap - 2d xsize-by-ysize array, snap shot matrix with size xsize
%                and ysize.
% 
% Output-only modal analysis toolbox v0.0
% Copyright by Hewenxuan Li, hewenxuan_li@uri.edu
% February 2022
snap = reshape(vec, xsize, ysize);
return