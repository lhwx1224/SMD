function vfill2curves(x, curve1, curve2, linecolor, fillcolor, facealpha)

if nargin > 3
    % Set the linecolors (upper and lower bound)
    LineColorUp = linecolor(1);
    LineColorDn = linecolor(2);
    FillColor = '#0071bd';
    FaceAlpha = 0.5;
    if nargin > 4
        FillColor = fillcolor;
        if nargin > 5
            FaceAlpha = facealpha;
        end
    end
else
    LineColorUp = 'k';
    LineColorDn = 'k';
    FillColor = '#0071bd';
    FaceAlpha = 0.5;
end
plot(x, curve1, LineColorUp, 'LineWidth', 1);
hold on;
plot(x, curve2, LineColorDn, 'LineWidth', 1);
x_contour = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x_contour, inBetween, hex2rgb(FillColor), 'FaceAlpha', FaceAlpha)
end

function [ rgb ] = hex2rgb(hex,range)
% Input checks:

assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.') 

if nargin==2
    assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end

% Tweak inputs if necessary: 

if iscell(hex)
    assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')
    
    % In case cell array elements are separated by a comma instead of a
    % semicolon, reshape hex:
    if isrow(hex)
        hex = hex'; 
    end
    
    % If input is cell, convert to matrix: 
    hex = cell2mat(hex);
end

if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end

if nargin == 1
    range = 1; 
end

% Convert from hex to rgb: 

switch range
    case 1
        rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;

    case {255,256}
        rgb = reshape(sscanf(hex.','%2x'),3,[]).';
    
    otherwise
        error('Range must be either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end

end