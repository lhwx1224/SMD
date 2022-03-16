function [delY, delTY] = TY2delY(TY, dx, dy, Operator)
% TY2delY(TY, dx, dy) differentiates a given spatiotemporal data tensor and
% return both the gradient-time tensor and the reshaped gradient-time
% field estimate.
%
% input: TY - tensor with dimension mm-by-nn-by-m, system response tensor of 
%        dx - sampling interval in x direction
%        dy - sampling interval in y direction
%        Operator - string, 'GradMag': magnitude of the gradient field
%                           'Laplacian': Laplacian operator (divergence field)
% output: delY - the spatially operated field measurement
%         delTY - the spatially operated field measurement tensor
%
% Copyright by Hewenxuan Li, hewenxuan_li@uri.edu
% Created on 2/15/2022 for Output-only Modal Analysis Toolbox v0.0
% Modified on 2/22/2022: added new operator "GradMag" the magnitude of the
%                        gradient field

if nargin < 4
    Operator = 'Laplacian';
end

% Inquire the size of the tensor
[mm, nn, m] = size(TY);

% Allocate memory for the gradient estimator
delY = zeros(m, mm*nn);
delTY = zeros(size(TY));

% FOR every temporal snapshot, DO:
for i = 1:m
    Wxy = TY(:,:,i);
    if isequal(lower(Operator), 'laplacian')
        [~, del2Wxy] = GenFiniteDiff2(Wxy, dx, dy);
        delY(i,:) = reshape(del2Wxy, 1, mm*nn);
        delTY(:,:,i) = del2Wxy;
    elseif isequal(lower(Operator), 'gradmag')
        [delWxy, ~] = GenFiniteDiff2(Wxy, dx, dy);
        delY(i,:) = reshape(delWxy, 1, mm*nn);
        delTY(:,:,i) = delWxy;
    else
        error('Check your input for the variable Operator; it should eigher be Laplacian or GradMag!')
    end
    
end