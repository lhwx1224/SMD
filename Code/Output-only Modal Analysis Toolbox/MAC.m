function [C, ModeIDIndex, SortIndex] = MAC(Phi_hat, Phi_a)
% MAC returns the modal assurance criterion matrix C and the sorted index
% for the analytical mode that is corresponding to the estimated mode
% shape. If no output variables specified, the matrix plot and a 3-d bar
% plot will be returned, indicating the quality of the modal
% identification.
% 
% Syntax: [C, ModeIDIndex, SortIndex] = MAC(Phi_hat, Phi_a)
%
% Input:    Phi_hat - estimated mode shapes
%             Phi_a - analytical/FEM mode shapes
%
% Output:         C - MAC matrix, correlation between estimated modes and
%                     the reference modes
%       ModeIDIndex - index indicating the closest reference mode shape to
%                     the currently considered esimated modes.
%         SortIndex - index indicating how the estimated modes should be
%                     sorted to follow the reference mode order.
%
%
% Copyright by Hewenxuan Li, 2021 
% References: [1] Modal Assurance Criterion - Pastor et. al., 2012

% Calculate the vector norm to each mode for both analytical and estimated.
norm_hat = vecnorm(Phi_hat);
norm_a = vecnorm(Phi_a);

% Make them squared
norm_hat = norm_hat.^2; % First term of the denominator
norm_a = norm_a.^2;     % Second term of the denominator

% Span the denominator matrix
D = norm_hat'*norm_a;

% Calculate the cross-covariance of the estimated modes and the analytical
% ones
N = Phi_hat'*Phi_a;

% Make it squared
N = N.^2;

% MAC Matrix
C = N./D;

% Index indicating the closest reference mode to the estimated modes
[~, ModeIDIndex] = max(C, [], 2);
% Index indicating the closest estimated mode to the reference modes
[~, SortIndex] = max(C, [], 1);

if nargout == 0
%     figure
    subplot(121)
    imagesc(C)
    colorbar
    % DEFINE THE VIRIDIS COLOR MAP
    color1 = '#FDE725';
    color2 = '#7AD151';
    color3 = '#22A884';
    color4 = '#2A788E';
    color5 = '#414487';
    color6 = '#440154';  
    map = customcolormap([0 0.2 0.4 0.6 0.8 1], {color1, color2, color3, color4, color5, color6});
    colormap(map)
    pbaspect([1 1 1])
    xlabel('Index of ${\Phi}$')
    ylabel('Index of $\hat{\Phi}$')
    subplot(122)
    h = bar3(C);
    % -------------- CODE FOR INTERPOLATED BAR COLOR ---------------------
    shading interp
    for i = 1:length(h)
        zdata = get(h(i),'Zdata');
        set(h(i),'Cdata',zdata)
        set(h,'EdgeColor','k')
    end
    set(h,'EdgeColor', 'none','LineWidth', 0.1);
% -------------------- CODE FOR SOLID BAR COLOR --------------------------
%     numBars = size(C,1);
%     numSets = size(C,2);
%     for i = 1:numSets
%         zdata = ones(6*numBars,4);
%         k = 1;
%         for j = 0:6:(6*numBars-6)
%             zdata(j+1:j+6,:) = C(k,i);
%             k = k+1;
%         end
%         set(h(i),'Cdata',zdata)
%     end
%     
% end
xlabel('Index of ${\Phi}$')
ylabel('Index of $\hat{\Phi}$')
zlabel('MAC')
pbaspect([1 1 1])
axis tight
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])
set(gcf, 'renderer', 'opengl')
end