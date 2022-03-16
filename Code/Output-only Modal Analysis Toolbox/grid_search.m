function param_matrix = grid_search(param_cell)
% PARAM_MATRIX = GRID_SEARCH(PARA_CELL) generates the matrix of parameters
% for parameter space exploration. The input is a cell array whose element
% contain a vector of unique values of parameter in one dimension. 
%
% syntax: param_matrix = grid_rearch(param_cell)
%
% input: param_cell, cell array whose element contains vector of uniqe
%        parameters in its dimension.
% output: param_matrix, matrix whose columns contains the location of the
%        first dimension of the hyperparameter space (hypercube).
%
% Hewenxuan Li, Dec 2020
x = param_cell{1};
y = param_cell{2};
% z = param_cell{3};
% q = param_cell{4};

[X,Y] = ndgrid(x,y);
xgrid = reshape(X,[],1);
ygrid = reshape(Y,[],1);
% zgrid = reshape(Z,[],1);
% qgrid = reshape(Q,[],1);
% param_matrix = [xgrid, ygrid, zgrid, qgrid];
param_matrix = [xgrid, ygrid];
end