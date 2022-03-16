%% SECTION 1: EIGENFUNCTIONS TO THE MEMBRANE VIBRATION PROBLEM (VISUALIZATION)
clear
a = 1;
b = 2;
rho = 1;
x = 0:0.01:a;
y = 0:0.01:b;

[Xgrids, Ygrids] = meshgrid(x,y);

m = 5;
n = 5;
W = cell(m,n);
for i = 1:m
    for j = 1:n
        W{i,j} = 2/sqrt(rho*a*b)*sin(i*pi*Xgrids/a).*sin(j*pi*Ygrids/b);
    end
end

% Put indices into cell arrays
param_cell{1} = 1:1:m;
param_cell{2} = 1:1:n;
% Grid search (all possible permutations of grids)
Mp = grid_search(param_cell);

clf
for k = 1:size(Mp,1)
    subplot(m,n,k)
    s = surf(W{Mp(k,1),Mp(k,2)});
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    pbaspect([b/a 1 1])
    axis tight
    title(['$W_{',num2str(Mp(k,1)),num2str(Mp(k,2)),'}$'],'Interpreter','latex')
    view(0,90)
end

%% SECTION 2: RITZ DISCRETIZATION
syms x y
rho0 = 1;
T0 = 1;
a0 = 1;
b0 = 1;
m = 1:1:2;
n = 1:1:2;
% Define the mass distribution
rho = 2 - (x - a0/2).^2 - (y - b0/2).^2;
% Define the eigenfunctions
W1 = sin(m*pi*x/a0);
W2 = sin(n*pi*y/b0);
W = 2/sqrt(rho0*a0*b0)*W1.'*W2;

% Put indices into cell arrays
param_cell{1} = 1:1:m(end);
param_cell{2} = 1:1:n(end);
% Grid search (all possible permutations of grids)
Mp = grid_search(param_cell);

% Discretize the eigenfunctions
M = zeros(m(end)*n(end),m(end)*n(end));
for k = 1:size(Mp,1)
    for l = 1:size(Mp,1)
        fun_m = matlabFunction(rho*W(Mp(k,1),Mp(k,2))*W(Mp(l,1),Mp(l,2)));
        M(k, l) = integral2(fun_m, 0, a0, 0, b0, "AbsTol", 1e-3);
    end
end

K = zeros(m(end)*n(end),m(end)*n(end));
for k = 1:size(Mp,1)
    for l = 1:size(Mp,1)
        g = gradient(W(Mp(k,1),Mp(k,2)),[x, y]).'*gradient(W(Mp(l,1),Mp(l,2)),[x, y]);
        fun_k = matlabFunction(T0*g);
        K(k, l) = integral2(fun_k, 0, a0, 0, b0, "AbsTol", 1e-3);
    end
end

%% SECTION 3: AMM RECTANGULAR MEMBRANE DEMO
mm = 3;
nn = 3;
dx = 0.01;
dy = 0.01;
a0 = 1;
b0 = 2;
rho0 = 1;
T0 = 1;
rho_xy_str = '2*rho - (x - a/2).^2 - (y - b/2).^2';
% rho_xy_str = 'rho*heaviside(x)*heaviside(y)';
T_xy_str = 'T*heaviside(x)*heaviside(y)';
F_xy_str = '1';
[M, K, U, Lambda, Phi, Phi_tilde, fw, omega_n] = ...
    AMM_uniform_membrane(mm, nn, dx, dy, a0, b0, rho0, T0, rho_xy_str, T_xy_str, F_xy_str);