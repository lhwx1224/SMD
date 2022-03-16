function Y = FirstLast_normalize(X)
% FIRSTLAST_NORMALISE(X) normalizes a vector or a matrix according to its
% last value and its first value in each of the vector (column-wise).
% syntax: Y = FirstLast_normalize(X)
%
% input: X - 1-D or 2-D array, input data which needs to be normalized
%
% output: Y - 1-D or 2-D array, same dimension to X, output normalized data
% 
% Hewenxuan Li Jan, 2021
if size(X,1) ~= 1 && size(X,2) ~= 1 
    Y = (X - X(1,:))./(X(end,:)-X(1,:));
else
    Y = (X - X(1))/(X(end)-X(1));
end
end