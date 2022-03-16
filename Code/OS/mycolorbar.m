function mycolormap = mycolorbar(MapName, Nscale)
% maps = mycolorbar(MapName, Nscale) renders the current graphics with customized
% colormap. When an output is designated, it outputs the colormap for
% future use.
% 
% Syntax: maps = mycolorbar(MapName, Nscale)
%
% Input: MapName - string, the name of the colormap, choose one from the
%                  list down below
%        Nscale - int, the total number of sample points in the color map,
%                 better if it is greater than 64.
% Output: maps - 2d array, RGB ratio triplet for the use of colormapping.
%
% MapNames: 'Viridis', 'BlackBody'
%
% Hewenxuan Li 2021
if nargin == 0
    MapName = 'Viridis';
    Nscale = [];
elseif nargin == 1
    Nscale = [];
end



switch lower(MapName)
    case 'viridis'
        % DEFINE THE VIRIDIS COLOR MAP
        color1 = '#FDE725';
        color2 = '#7AD151';
        color3 = '#22A884';
        color4 = '#2A788E';
        color5 = '#414487';
        color6 = '#440154';
        if isempty(Nscale)
            mycolormap = customcolormap([0 0.2 0.4 0.6 0.8 1], {color1, color2, color3, color4, color5, color6});
            if nargout == 0
                colormap(mycolormap);
            end
        else
            mycolormap = customcolormap([0 0.2 0.4 0.6 0.8 1], {color1, color2, color3, color4, color5, color6}, Nscale);
            if nargout == 0
                colormap(mycolormap);
            end
        end

    case 'blackbody'
        if isempty(Nscale)
            mycolormap = customcolormap(linspace(0,1,15), {'#000000', '#240f09', '#3e1611',...
                '#5a1b16', '#771e1a', '#96211e', '#b42622', '#c5411c', '#d65813', '#e47007',...
                '#e78d12', '#e9c327', '#e7de32', '#f6f090', '#ffffff'});
            if nargout == 0
                colormap(mycolormap)
            end
        else
            mycolormap = customcolormap(linspace(0,1,15), {'#000000', '#240f09', '#3e1611',...
                '#5a1b16', '#771e1a', '#96211e', '#b42622', '#c5411c', '#d65813', '#e47007',...
                '#e78d12', '#e9c327', '#e7de32', '#f6f090', '#ffffff'}, Nscale);
            if nargout == 0
                colormap(mycolormap)
            end
        end
end