function [delX, del2X]= GenFiniteDiff2(X, dx, dy)
% GenFiniteDifference calculates the derivatives up to the second order,
% given a data matrix X, the differentiation step dx, and the options which
% designate the methodology to be used.
% 
% Syntax: [dX, ddX] = GenFiniteDiff(X, dx, dy)
%
% Input: X - data matrix 
%        dx - step between two samples of data for differentiation
%        dy
%
% Output: X - trimmed data matrix
%         dX - 1st order derivative of X
%         ddX - 2nd order derivative of X
%
% Copyright by Hewenxuan Li, June 2021, hewenxuan_li@uri.edu
% 
% Generated for HOD project

% Check if input is a matrix
if size(X,1) == 1 || size(X,2) == 1
    error('Input must be a matrix!')
end

% memory allocation
[m, n] = size(X);
delX = zeros(m, n);
delXx = zeros(m, n);
delXy = zeros(m, n);
del2X = zeros(m, n);


% Calculate the first-order derivative (discrete Del operator)
% For inner nodes DO:
for i = 2:m-1
    for j = 2:n-1
        delXx(i,j) = (X(i + 1, j) - X(i - 1, j))/2/dx; % Grad. in x
        delXy(i,j) = (X(i, j + 1) - X(i, j - 1))/2/dy; % Grad. in y
    end
end
% For boundaries DO: (Dirichlet - forward/backward approximation)
% Note: developed for general field, if regular field under consideration,
% vectorial construction will be faster
for i = 1
    for j = 1:n
        if j == 1
            delXx(i,j) = (X(i + 1, j) - X(i, j))/dx; % (1,1) corner in x
            delXy(i,j) = (X(i, j + 1) - X(i, j))/dy; % (1,1) corner in y
        elseif j == n
            delXx(i,j) = (X(i + 1, j) - X(i, j))/dx; % (1,n) corner in x
            delXy(i,j) = (X(i, j) - X(i, j - 1))/dy; % (1,n) corner in y
        else
            delXx(i,j) = (X(i + 1, j) - X(i, j))/dx; % (1,j) internal points in x
            delXy(i,j) = (X(i, j + 1) - X(i, j - 1))/2/dy; % (1,j) internal points in y
        end
    end
end

for i = m
    for j = 1:n
        if j == 1
            delXy(i,j) = (X(i, j) - X(i - 1, j))/dx; % (m,1) corner in x
            delXy(i,j) = (X(i, j + 1) - X(i, j))/dy; % (m,1) corner in y
        elseif j == n
            delXx(i,j) = (X(i, j) - X(i - 1, j))/dx; % (m,n) corner in x
            delXy(i,j) = (X(i, j) - X(i, j - 1))/dy; % (m,n) corner in y
        else
            delXx(i,j) = (X(i, j) - X(i - 1, j))/dx; % (m,j) internal points in x
            delXy(i,j) = (X(i, j + 1) - X(i, j - 1))/2/dy;  % (m,j) internal points in y
        end
    end
end

% For the first column center diff in x forward diff in y
for j = 1
    for i = 2:m-1
        delXx(i,j) = (X(i+1, j) - X(i-1, j))/2/dx;
        delXy(i,j) = (X(i, j + 1) - X(i, j))/dy;
    end
end

% For the last column center diff in x backward diff in y
for j = n
    for i = 2:m-1
        delXx(i,j) = (X(i+1, j) - X(i-1, j))/2/dx;
        delXy(i,j) = (X(i, j) - X(i, j - 1))/dy;
    end
end

delX = sqrt(delXx.^2 + delXy.^2);

% ------------------------------------------------------------------------
% THE SECOND ORDER DERIVATIVE IS INDEPENDENT OF THE PREVIOUS STEPS
% ------------------------------------------------------------------------

% Calculate the second-order derivative for interior points 
% (discrete Laplace Operator)
for i = 2:m-1
    for j = 2:n-1
        del2X(i,j) = (X(i+1, j) - 2*X(i, j) + X(i-1, j))/dx^2 + (X(i, j + 1) - 2*X(i, j) + X(i, j - 1))/dy^2;
    end
end

% Calculate the Laplacian at boundary points (x-direction)
for i = 1
    for j = 2:n-1
        del2X(i,j) = (X(i, j) - 2*X(i + 1, j) + X(i + 2, j))/dx^2 + (X(i, j + 1) - 2*X(i, j) + X(i, j - 1))/dy^2;
    end
end

for i = m
    for j = 2:n-1
        del2X(i,j) = (X(i - 2, j) - 2*X(i - 1, j) + X(i, j))/dx^2 + (X(i, j + 1) - 2*X(i, j) + X(i, j - 1))/dy^2;
    end
end

% Calculate the Laplacian at the boundary points (y-direction)
for i = 2:n-1
    for j = 1
        del2X(i,j) = (X(i + 1, j) - 2*X(i, j) + X(i - 1, j))/dx^2 + (X(i, j + 2) - 2*X(i, j + 1) + X(i, j))/dy^2;
    end
end

for i = 2:n-1
    for j = n
        del2X(i,j) = (X(i + 1, j) - 2*X(i, j) + X(i - 1, j))/dx^2 + (X(i, j) - 2*X(i, j - 1) + X(i, j - 2))/dy^2;
    end
end