function [pvec, pmark] = plotvecs(v1, v2, PlotColor, LineStyle, Marker)
% plotvecs.m plots vectors using two data points that are associated to the
% tail and the head of the vector.
%
% syntax: [pvec, pmark] = plotvecs(v1, v2, PlotColor, LineStyle, Marker)
%
% input: v1 - vector of the tail (origin) of the vector being plotted
%        v2 - vector of the head of the vector being plotted
%        PlotColor - color of the vector being plotted, e.g., 'r'
%        LineStyle - line style of the vector being plotted, e.g., '-'
%        Marker - marker of the vector being plotted that indicates the
%                 head of the vector.
%
% Hewenxuan Li 2021
% Output-only Modal Analysis Toolbox v0.0

[m1, n1] = size(v1);
[m2, n2] = size(v2);

% if the input is a row vector, turn it into a column vector
if m1 < n1
    v1 = v1';
end
if m2 < n2
    v2 = v2';
end

% Check the vector size to see if it is a 2d or 3d one
[m1, n1] = size(v1);
[m2, n2] = size(v2);

% Check if the vectors are vectors and if they have identical dimensions
if n1 ~= 1 || n2 ~= 1
    error('Input must be a vector!')
end

if m1 ~= m2
    error('Vectors must have the same dimension!')
end

% Check if the vector is of dimension 2 or 3
if m1 ~= 2 && m1~= 3
    error('The first vector has incompatible size to be plotted!')
end

if m2 ~= 2 && m2~= 3
    error('The second vector has incompatible size to be plotted!')
end

% column partition of V is the two vectors
V = [v1 v2];

% extract vectors to be plotted
if m1 == 2
    v1p = V(1,:);
    v2p = V(2,:);
    if nargin < 3
        pvec = plot(v1p, v2p);
        pmark = [];
    else
        pvec = plot(v1p, v2p, "Color",PlotColor,"LineStyle",LineStyle);
        if nargin == 5
            hold on
            pmark = plot(v1p(end), v2p(end), "Color",PlotColor,"Marker",Marker);
        else
            pmark = [];
        end
    end
elseif m1 == 3
    v1p = V(1,:);
    v2p = V(2,:);
    v3p = V(3,:);
    if nargin < 3
        pvec = plot3(v1p, v2p, v3p);
        pmark = [];
    else
        pvec = plot3(v1p, v2p, v3p, "Color",PlotColor,"LineStyle",LineStyle);
        if nargin == 5
            hold on
            pmark = plot3(v1p(end), v2p(end), v3p(end), "Color",PlotColor,"Marker",Marker);
        else
            pmark = [];
        end
    end
end


