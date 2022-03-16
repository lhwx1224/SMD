function [som, sov_tsmd, spm_tsmd, sc_tsmd, S1_tsmd, S2_tsmd, V1_tsmd, V2_tsmd] = tsmd2s(Y, delY, r, MeanSubtraction)
% ------------------------------------------------------------------------
%        Transformed/Truncated Smooth Mode Decomposition (Spatial)
%      --------------------------------------------------------------
% 1. Temporal dimensionality reduction by obtaining an orthogonal
%    projector, POC of [X; delX], and truncated it to some degree
% 2. Data projection onto the common POM, Ut
% 3. SMD using the projected data sets [Xt; delXt]
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
%                          Yt = Ut'*Y
%                       delYt = Ut'delY
%
%           Conduct GSVD to the temporally tranformed matrix pair:
%
%                          Yt' = Usct*S1t*V1t'
%                       delYt' = Usct*S2t*V2t'   <Truncated Smooth Coords>
%
%            Transform Usct back to the original dimensionality, Usc:
% 
%                          Usc = Ut*Usct;        <Smooth Coordinates>
% ------------------------------------------------------------------------
%
% syntax: [sc, sov_tsmd, spm_tsmd, som_tsmd, S1_tsmd, S2_tsmd, V1_tsmd, V2_tsmd] = tsmd2(Y, delY, r, MeanSubtraction)
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
% Created on 2/22/2022, Output Only Modal Analysis Toolbox v0.0
% ------------------------------------------------------------------------

% ---------------------- CHECK INPUT VARIABLES --------------------------
[m, n] = size(Y); % Check the size of the input field matrix
if m > n
    error('TSMD expects a matrix with more columns than its rows; i.e., m < n with [m, n] = size(Y)!')
end

% Check the input variables and set the variables to default values
if nargin < 4
    MeanSubtraction = 'False';
    warning('Mean of the data will not be subtracted from the data during SOD!')
elseif nargin < 3
    MeanSubtraction = 'False';
    warning('Mean of the data will not be subtracted from the data during SOD!')
    r = m;
    warning('No truncation will be applied! r = m!')
end

% ---------------- 1.DIRECT POD TO [X, delX] ----------------------------

[~, St, Vt] = svd([Y;delY]); % <Economy not in use to resolve all possible subspaces>

% --------- 2.TSVD-BASED TEMPORAL DIMENSIONALITY REDUCTION ---------------

Q = Vt(:, 1:r);  % Orthogonal basis by truncating POCs 
Yr = Y*Q;       % Temporally-reduced measurement matrix
delYr = delY*Q; % Temporally-reduced gradient matrix

% ---------- 3. CONDUCT SMOOTH MODE DECOMPOSITION ------------------------

[som_tsmd, sov_tsmd, spm_tsmd, sc_tsmd, S1_tsmd, S2_tsmd, V1_tsmd, V2_tsmd] = sod(Yr', delYr', ' ', MeanSubtraction);

% ---------- 4. TRANSFORM BACK TO THE ORIGINAL DIMENSION -----------------

som = Q*som_tsmd;

% ------------- VISUALIZE THE POVS OF THE TRANSFORMATION -----------------

figure(64),clf
plot(diag(St), '--o')
pbaspect([2 1 1])
xlabel('Index of Propoer Orthogonal Subspaces')
ylabel('POV($[Y, \mathcal{D}_s(Y)]$)')
set(gca, 'Yscale', 'log')
title('Proper Orthogonal Values of $[Y, \mathcal{D}_s(Y)]$')
grid on

return