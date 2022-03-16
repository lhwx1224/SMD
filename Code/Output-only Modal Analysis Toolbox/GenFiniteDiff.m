function [dX, ddX]= GenFiniteDiff(X, dx, options)
% GenFiniteDifference calculates the derivatives up to the second order,
% given a data matrix X, the differentiation step dx, and the options which
% designate the methodology to be used.
% 
% Syntax: [dX, ddX] = GenFiniteDiff(X, dx, options)
%
% Input: X - data matrix 
%        dx - step between two samples of data for differentiation
%        options - string, 'f1' - first order forward difference
%                          'f2' - second order forward difference
%                          'c2' - second order center difference
%                          'c4' - fourth order center difference
% Output: X - trimmed data matrix
%         dX - 1st order derivative of X
%         ddX - 2nd order derivative of X
%
% Copyright by Hewenxuan Li, June 2021, hewenxuan_li@uri.edu
% 
% Generated for HOD project

% memory allocation
[m, n] = size(X);
dX = zeros(m, n);
ddX = zeros(m, n);
% Case switch
switch lower(options)
    case 'f1'
        % O(dx) forward finite difference
        for i = 1:m-1
            dX(i,:) = X(i+1,:) - X(i, :);
        end
        dX(m, :) = [];
        dX = dX/dx;
        ddX = [];
        X = X(1:end-1,:); % Trim the original data
        
    case 'f2'
        % O(dx^2) forward finite difference
        for i = 1:m-2
            dX(i,:) = (- 3*X(i, :) + 4*X(i+1,:) - X(i+2,:));
        end
        dX = dX/2/dx;
        dX(m-2:m,:) = []; % Trim the derivative data
        X = X(1:end-3,:); % Trim the original data

        for i = 1:m-3
            ddX(i,:) = (2*X(i, :) - 5*X(i+1,:) + 4*X(i+2,:) - X(i+3,:));
        end
        ddX = ddX/(dx^3);
        ddX(m-2:m,:) = [];
        
    case 'c2'
        % O(dx^2) centered finite difference
        % TESTING SPLINEDIFF
%         dX = splinediff(X, 1/dx);
        for i = 2:m-1
            dX(i,:) = 1/2*(X(i+1,:) - X(i-1,:));
        end
%         dX([1, m],:) = [];
        
        % Consider Edge Effects
        % LEGACY CODE ------------------------------
%         for i = 1
%             X0 = X(2,:) - 2*dX(2,:);
%             dX1 = (X(2,:) - X0)/2;
%             dX(i, :) = dX1;
%         end
% 
%         for i = m
%             Xmp1 = X(m-1,:) + 2*dX(m-1,:);
%             dXm = (Xmp1 - X(m-1,:))/2;
%             dX(m,:) = dXm;
%         end
        % ------------------------------------------
        
        % Consider Edge Effects (using forward or backward differentiation)
        for i = 1
            dX(i,:) = X(2,:) - X(i,:);
        end

        for i = m
            dX(i,:) = X(m,:) - X(m - 1,:);
        end

        dX = dX/dx;

        % O(dx^2) centered finite difference 2nd order
        for i = 2:m-1
            ddX(i,:) = X(i+1,:) - 2*X(i,:) + X(i-1,:);
        end
        ddX([1, m],:) = [];
        X = X(2:m-1,:); % Trim the original data
    case 'c4'
        % O(dx^4) centered difference approximations
        for i = 3:m-2
            dX(i,:) = -X(i+2, :) + 8*X(i+1,:) - 8*X(i-1,:) + X(i-2,:);
        end
        dX([1,2,m-1,m],:) = [];
        dX = dX/12/dx;
        for i = 3:m-2
            ddX(i,:) = -X(i+2) + 16*X(i+1) -30*X(i,:) + 16*X(i-1,:) + X(i-2,:);
        end
        ddX([1,2,m-1,m],:) = [];
        ddX = ddX/12/(dx^2);
        X = X(3:m-2,:); % Trim the original data
end

