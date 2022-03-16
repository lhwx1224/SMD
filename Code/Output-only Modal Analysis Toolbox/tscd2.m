function [sm, sov_tscd, spm_tscd, soc_tscd, S1_tscd, S2_tscd, U1_tscd, U2_tscd] = tscd2(Y, DY, delY, r, MeanSubtraction)
% ------------------------------------------------------------------------
%        Transformed/Truncated Smooth Coordinate Decomposition (Spatial)
%      --------------------------------------------------------------
% 1. Temporal dimensionality reduction by obtaining an orthogonal
%    projector, POC of [X; delX], and truncated it to some degree
% 2. Data projection onto the common POM, Ut
% 3. SCD (SOD) using the projected data sets [Xt, delXt]
% 4. Unfold the modes to its original dimensionality 
%
%   TSMD:  given Y and delY, conduct SVD([Y; delY]) to find the common
%          orthognal coordinates U:
%               
%                          Y = U1*Sigma*V'
%                       delY = U2*Sigma*V'
%
%          Obtain transformation matrix for temporal reduction:
%
%                         Vt = V(:, 1:r);
%          
%          Temporal dimensionality reduction by using Ut as the input space
%          to sample the rows of the matrix Y and delY:
%
%                          Yt = Y*Vt
%                       delYt = delY*Vt
%
%           Conduct GSVD to the temporally tranformed matrix pair:
%
%                          Yt' = U1t*S1t*Vsmt'
%                       delYt' = U1t*S2t*Vsmt'  <Truncated Smooth Modes>
%
%            Transform Usct back to the original dimensionality, Usc:
% 
%                          Vsm = Vt*Vsmt;    <Smooth Modes>
% ------------------------------------------------------------------------
%
% syntax: [sm, sov_tscd, spm_tscd, soc_tscd, S1_tscd, S2_tscd, U1_tscd, U2_tscd] = tsmd(Y, DY, delY, r, MeanSubtraction)
%
% input: Y - 2d array, Flattened field measurement data
%        delY - 2d array, Flattended spatially operated measurement data
%        r - integer, degree of truncation (default is # of rows of Y)
%        MeanSubtraction - string, subtract the mean from the data set
%                          (1) - 'True', subtract mean from data
%                          (2) - 'False', no mean subtraction
%
% output: sc - 2d array of size m-by-r, smooth coordinate from TSMD
%         sov_tsmd - 1d array of size r-by-1, smooth orthogonal value in
%                    reduced space.
%         spm_tsmd - 2d array of size r-by-r, smooth projective modes in
%                    reduced space.
%         som_tsmd - 2d array of size r-by-r, smooth orthogonal modes in
%                    reduced space.
%         S1_tsmd - 1d array of size r-by-1, diagnonal element of S1 in
%                   GSVD.
%         S2_tsmd - 1d array of size r-by-1, diagnonal element of S2 in
%                   GSVD.
%         V1_tsmd - 2d array of size r-by-r, diagnonal element of V1 in
%                   GSVD.
%         V1_tsmd - 2d array of size r-by-r, diagnonal element of V2 in
%                   GSVD.
%
%
% Copyright by Hewenxuan Li, hewenxuan_li@uri.edu
% Created on 2/20/2022, Output Only Modal Analysis Toolbox v0.0
% ------------------------------------------------------------------------

% ---------------------- CHECK INPUT VARIABLES --------------------------
[~, n] = size(Y); % Check the size of the input field matrix

% Check the input variables and set the variables to default values
if nargin < 5
    MeanSubtraction = 'False';
    warning('Mean of the data will not be subtracted from the data during SOD!')
elseif nargin < 4
    MeanSubtraction = 'False';
    warning('Mean of the data will not be subtracted from the data during SOD!')
    r = n;
    warning('No truncation will be applied! r = n!')
elseif nargin < 3
    error('Input should at least involve Y, DY, and delY!')
end

% ---------------- 1.DIRECT POD TO [X; delX] ----------------------------

[~, St, Vt] = svd([Y; delY], 'econ');

% --------- 2.TSVD-BASED SPATIAL DIMENSIONALITY REDUCTION ---------------

Q = Vt(:, 1:r);  % Orthogonal basis by truncating POCs 
Yt = Y*Q;       % Temporally-reduced measurement matrix
DYt = DY*Q; % Temporally-reduced gradient matrix

% ---------- 3. CONDUCT SMOOTH MODE DECOMPOSITION ------------------------

[soc_tscd, sov_tscd, spm_tscd, sm_tscd, S1_tscd, S2_tscd, U1_tscd, U2_tscd] = sod(Yt, DYt, ' ', MeanSubtraction);

% ---------- 4. TRANSFORM BACK TO THE ORIGINAL DIMENSION -----------------

sm = Q*sm_tscd;

% ------------- VISUALIZE THE POVS OF THE TRANSFORMATION -----------------

figure(64),clf
plot(diag(St), '--o')
pbaspect([2 1 1])
xlabel('Index of Propoer Orthogonal Subspaces')
ylabel('POV($[Y, \mathcal{D}_s(Y)]$)')
set(gca, 'Yscale', 'log')
title('Proper Orthogonal Values of $[Y; \mathcal{D}_s(Y)]$')
grid on

return