function [soc, sov, spm, som, C, S, U, V] = sod(Y,DY,DiffOpt,MeanSub)
% Smooth Orthogonal Decomposition through generalized singular value
% decomposition (GSVD).
% sod takes in a 2d array Y whose columns are acquired data, each column
% stores a univariate time series (has to be time series!). If only one
% data set is provided, the program will automatically calculate the
% temporal differentiated data set DY through designated differentiation
% method.
% 
% syntax: [soc, sov, spm, som, C, S, U, V] = sod(Y,DY,DiffOpt,MeanSub)
%
% input: Y  - 2d-array, original input data matrix
%        DY - 2d-array, the smooth-operated data matrix via smoothness
%             operator D(), the default operator is the temporal
%             differentiation operator that approximates \dot{Y}.
%        DiffOpt - string, name of the differential operator
%                  ' ' - when a space is used, the finite difference
%                        operator using MATLAB diff function is applied
%                  'spectral' - spectral diffentiation 
%                  'spline' - spline-fit-based differentiation
%        MeanSub - string, condition that stipulate if the mean value of
%                  the data is subtracted from the raw data
%                  'True' - subtract the mean from the data 
%                  'False' - keep the raw data as is
%
% output: soc - smooth orthogonal coordinates (SOCs), 2d-array, obtained by
%               projecting the data set Y onto the smooth projective modes
%               (SPMs) using soc = Y*spm
%         sov - smooth orthogonal values (SOVs), 1d-array that stores all
%               the smooth orthogonal values calculated as the ratios of
%               the singular values of the form sov_i = c_i/s_i
%         spm - smooth projective modes (SPMs), 2d-array, obtained by
%               finding the contravariant basis of the resolved smooth
%               modes (SMs), spm = inv(som')
%         sm  - smooth modes, 2d-array, obtained by simultaneously
%               decompose Y and DY via GSVD of the form: 
%
%                         Y = U*C*X' and DY = V*S*X'
%
%               som = fliplr(X);
%           C - singular values matrix to the first decomposition
%           S - singular values matrix to the second decomposition
%           U - left singular matrix to the first decomposition, 2d-array,
%               that stores normalized smooth orthogonal coordinates
%           V - left singular matrix to the second decomposition, 2d-array,
%               that stores normalized velocity coordinates
%
% Created by Hewenxaun Li, May 2021
% -------------------------------------------------------------------------


if nargin<3
    MeanSub = 'True';
elseif nargin<2
    MeanSub = 'True';
    DiffOpt = ' ';
end
%% initialize parameters and center the data matrix
[n, ~] = size(Y); % n--temporal dimension; m--spatial dimension
% Mean Subtraction
if isequal(lower(MeanSub),'true')
    Y = Y - ones(n,1)*mean(Y); % center data columns (make them zero mean)
end

if nargin == 3
    if isequal(DiffOpt,'spline')
        DY = splinediff(Y);
    elseif isequal(DiffOpt,'spectral')
        DY = spectraldiff(Y);
    else
        DY = diff(Y);
    end
end
    
%% perform the generalized sod
[U,V,X,C,S] = gsvd(Y,DY,0);

if size(C,1) > size(S,1)
    sov = flipud(diag(C(1:size(S,1),1:size(S,1)))./diag(S(1:size(S,1),1:size(S,1))));
else
    sov = flipud(diag(C)./diag(S));
end
spm = fliplr(inv(X'));
som = fliplr(X);
soc = Y*spm;
C = flip(diag(C));
S = flip(diag(S));
U = fliplr(U);
V = fliplr(V);
end