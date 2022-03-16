function TY = Y2TY(Y, mm, nn)
% Y2TY transforms a data matrix whose colummns are vectors reshaped from
% spatial snapshots. The jth column of Y corresponds to the jth spatial
% snapshot of the spatiotemporal tensor.
%
% syntax: [TY, Nx, Ny, Nt] = Y2TY(Y, mm, nn)
%
% input: Y - m-by-n 2d array whose columns are vectors reshaped from
%            spatial snapshots whose x has mm elements and whose y has nn
%            elements
%        mm - # of rows in the orignial snapshot
%        nn - # of columns in the original snapshot
%
% output: TY - reshaped and temporally organized tensor whose first two
%              dimensions correspond to the snapshot of the original field
%
% Copyright by Hewenxuan Li, hewenxuan_li@uri.edu
%
% Created on 2/15/2022

% Inquire the size of the data matrix
[m, ~] = size(Y);

% Tensor with mm-by-nn spatial snapshot and size(Y,2) depth in time
TY = zeros(mm, nn, m); 

% FOR every time step, DO:
for i = 1:m % Assuming spatial dimension is distributed column-wise 
    TY(:,:,i) = reshape(Y(i,:), mm, nn);
end

return

