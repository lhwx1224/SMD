function [M, K, U, Lambda, Phi, Phi_tilde, fw, omega_n] = AMM_uniform_FF_membrane(mm, nn, dx, dy, a0, b0, rho0, T0, rho_xy_str, T_xy_str, F_xy_str)
% AMM_uniform_membrane() discretizes a rectangular membrane with fixed
% boundary condition into mass, stiffness matrices, and a weighting vector
% to the forcing function in the configuration space using the eigenfunctions 
%
% syntax: [M, K, U, Lambda, Phi, Phi_tilde, fw, omega_n] = AMM_uniform_membrane(mm, nn, dx, dy, a0, b0, rho0, T0, rho_xy_str, T_xy_str, F_xy_str)
%
% input: mm - number of eigenfunctions in the x direction
%        nn - number of eigenfunctinos in the y direction
%        dx - spatial resolution for the x eigenfunction
%        dy - spatial resolution for the y eigenfunction
%        a0 - edge legnth of the x direction
%        b0 - edge length of the y direction
%        rho0 - the mass density constant (default = 1)
%        T0 - the constant part of the tension (default = 1)
%        rho_xy_str - mass density function string description
%        T_xy_str - tension function string description
%
% output: M - discretized mass matrix
%         K - discretized stiffness matrix
%         U - modal transformation matrix for undamped system (that transforms the assumed modes to the ``true'' modes)
%         Lambda - modal dynamics matrix for undamped system (that embed the natural frequencies of the ``true'' system)
%         Phi - cell, transformed ``true'' modes for undamped system
%         Phi_tilde - cell, assumed modes
%         fw - force weight
%         omega_n - natural frequencies of the ``true'' modal system
%
% Example: mm = 2;
%          nn = 2;
%          dx = 0.01;
%          dy = 0.01;
%          a0 = 1;
%          b0 = 2;
%          rho0 = 1;
%          T0 = 1;
%          rho_xy_str = '2*rho - (x - a/2).^2 - (y - b/2).^2';
%          rho_xy_str = 'rho*heaviside(x)*heaviside(y)';
%          T_xy_str = 'T*heaviside(x)*heaviside(y)';
%          F_xy_str = '1';
%          [M, K, U, Lambda, Phi, Phi_tilde, fw, omega_n] = ...
%          AMM_uniform_membrane(mm, nn, dx, dy, a0, b0, rho0, T0, rho_xy_str, T_xy_str, F_xy_str);
%
%          *OR RUN SECTION 3 in Rectangular_Membrane_Demo.m
%
% Copyright by Hewenxuan Li, hewenxuan_li@uri.edu
% January 2022
%-------------------------------------------------------------------------
%% UNIFORM EULER BERNOULLI BEAM AND ITS EIGENFUNCTIONS   
disp("================= FIXED-FIXED MEMBRANE AMM DISCRETIZATION STARTED! =========================")
addpath('L:\My Drive\Graduate study\Research\Projects\OS')
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
warning('off')
% Check if the parpool is running; if not, execute:
if isempty(gcp('nocreate'))
    clusterinfo = parcluster('local');
    Nworkers = clusterinfo.NumWorkers;
    parpool('local', Nworkers)
end
tic
syms x y a b rho T
global x y
m = 1:1:mm;
n = 1:1:nn;
nm = mm*nn; % Total number of modes considered
fprintf('Number of eigenfunctions in x direction: %d.\n', mm)
fprintf('Number of eigenfunctions in y direction: %d.\n', nn)
fprintf('Total number of assumed modes considered: %d.\n', nm)

% Define the mass distribution
rho_xy = str2sym(rho_xy_str);
T_xy = str2sym(T_xy_str);
F_xy = str2sym(F_xy_str);

rho_xy = subs(rho_xy, [a b rho], [a0 b0 rho0]);
T_xy = subs(T_xy, [a b T], [a0 b0 T0]);
F_xy = subs(F_xy, [a b], [a0 b0]);

% Define the normalized eigenfunctions
W1 = sin(m*pi*x/a0);
W2 = sin(n*pi*y/b0);
W = 2/sqrt(rho0*a0*b0)*W1.'*W2;

figure(1), clf
xs = 0:dx:a0;
ys = 0:dy:b0;
[Xgrids, Ygrids] = meshgrid(xs,ys);
subplot(121)
fsf1 = fsurf(matlabFunction(rho_xy),[0 a0 0 b0]);
pbaspect([a0/b0 1 a0/2])
fsf1.EdgeColor = 'none';
title('Mass density distribution'),xlabel('$x$'),ylabel('$y$'),zlabel('$\rho(x,y)$')
mycolorbar();

subplot(122)
fsf2 = fsurf(matlabFunction(T_xy),[0 a0 0 b0]);
pbaspect([a0/b0 1 a0/2])
fsf2.EdgeColor = 'none';
title('Tension distribution'),xlabel('$x$'),ylabel('$y$'),zlabel('$T(x,y)$')
mycolorbar();

sgtitle('Mass Distribution and Tension Distribution of the Membrane')
%% ---------------- Eigenvalue Problem --------------------------------
% Put indices into cell arrays
param_cell{1} = 1:1:m(end);
param_cell{2} = 1:1:n(end);
% Grid search (all possible permutations of grids)
Mp = grid_search(param_cell);

% Discretize the eigenfunctions
M = zeros(m(end)*n(end),m(end)*n(end));
K = zeros(m(end)*n(end),m(end)*n(end));

% tic
% for k = 1:size(Mp,1)
%     for l = 1:size(Mp,1)
%         fun_m = matlabFunction(rho_xy*W(Mp(k,1),Mp(k,2))*W(Mp(l,1),Mp(l,2)));
%         M(k, l) = integral2(fun_m, 0, a0, 0, b0, "AbsTol", 1e-3);
%     end
% end
% 
% for k = 1:size(Mp,1)
%     parfor l = 1:size(Mp,1)
%         g = gradient(W(Mp(k,1),Mp(k,2)),[x, y]).'*gradient(W(Mp(l,1),Mp(l,2)),[x, y]);
%         fun_k = matlabFunction(T_xy*g);
%         K(k, l) = integral2(fun_k, 0, a0, 0, b0, "AbsTol", 1e-3);
%     end
% end
% t1 = toc;
% ------------------------------------------------------------------------
% PARALELLIZE THE INTEGRATION - SPEED UP: 38.4% (intel i9-9900k 4.7GHz)
% ------------------------------------------------------------------------
disp('------- Paralallizing the discretization for mass matrix - M ------')
for k = 1:size(Mp,1)
    parfor l = 1:size(Mp,1)
        fun_m = matlabFunction(rho_xy*W(Mp(k,1),Mp(k,2))*W(Mp(l,1),Mp(l,2)));
        M(k, l) = integral2(fun_m, 0, a0, 0, b0, "AbsTol", 1e-3);
    end
end
disp('--------------- M-discretization Completed! -----------------------')

disp('---- Paralallizing the discretization for stiffness matrix - K ----')
for k = 1:size(Mp,1)
    parfor l = 1:size(Mp,1)
        g = gradient(W(Mp(k,1),Mp(k,2)),[x, y]).'*gradient(W(Mp(l,1),Mp(l,2)),[x, y]);
        fun_k = matlabFunction(T_xy*g);
        K(k, l) = integral2(fun_k, 0, a0, 0, b0, "AbsTol", 1e-3);
    end
end
disp('--------------- K-discretization Completed! -----------------------')

% Calculate the weights to the forcing function
% Define a force variable f(x,t) = F(x)G(t) and we only worries about the
% spatial part F(x) since the integration will be conducted over space.
syms f
fw = zeros(nm, 1);
% Fx = dirac(x - L0); % <ANALYTICAL>
disp('--------------- Forcing Weight Paralillizing ----------------------')
parfor k = 1:nm
    fun_f = matlabFunction(F_xy*W(Mp(k,1),Mp(k,2)));
    fw(k) = integral2(fun_f, 0, a0, 0, b0, "AbsTol", 1e-3); % <ANALYTICAL>
%     progress_bar(k,nm,'Calculating Forcing Weights')
end
disp('----------------- Forcing Weight Completed! -----------------------')

% Solve for the system dynamics (Symmetric Eigenvalue Problem - EVP)
[evec, eval] = eig(inv(M)*K);
[eval_s, I] = sort(diag(eval));
eval = diag(eval_s);
evec = evec(:, I);
omega_n = sqrt(diag(eval));

%% Assign output parameters
U = evec;
Lambda = eval;
Phi_tilde = cell(size(Mp,1),1);
Phi = cell(size(Mp,1),1);
Phi_tilde_mat = zeros(length(xs)*length(ys), size(Mp,1));
parfor i = 1:nm
    % Assign the symbolic function to a matlab function handle
    fun_W1 = matlabFunction(W1(Mp(i,1)));
    fun_W2 = matlabFunction(W2(Mp(i,2)));
    % Evaluate the x-wise  and y-wise functions separately
    W_x = fun_W1(xs);
    W_y = fun_W2(ys);
    Phi_tilde{i} = W_x'*W_y;
    Phi_tilde_mat(:,i) = reshape(Phi_tilde{i},[],1);
    % Evaluate the function using the xs and ys value pairs
    % Do loop through the column-wise dimension
%     input_xy = [ones(length(ys),1)*xs(i), ys'];
%     Phi_tilde(:,i) = fun_W(input_xy);
end

Phi_mat = Phi_tilde_mat*U;

parfor i = 1:nm
    Phi{i} = reshape(Phi_mat(:,i),length(xs),length(ys));
end

time_taken = toc;

plotlim = 5; % plot maximally a plotlim-by-plotlim figures

figure(2),clf
for k = 1:min(mm*nn, plotlim^2)
    subplot(min(mm,plotlim),min(nn,plotlim),k)
    img = imagesc(Phi_tilde{k});
    pbaspect([b0/a0 1 1])
    mycolorbar();
    title(['$\widetilde{\Phi}_{',num2str(Mp(k,1)),num2str(Mp(k,2)),'}$'],'Interpreter','latex')
    xlabel('$x$')
    ylabel('$y$')
end
sgtitle(['Assumed Modes: $\widetilde{\Phi}$'])


figure(3),clf
for k = 1:min(mm*nn, plotlim^2)
    subplot(min(mm,plotlim),min(nn,plotlim),k)
    img = imagesc(Phi{k});
    pbaspect([b0/a0 1 1])
    mycolorbar();
    title(['${\Phi}_{',num2str(Mp(k,1)),num2str(Mp(k,2)),'}$'],'Interpreter','latex')
    xlabel('$x$')
    ylabel('$y$')
end
sgtitle('Approximated Modes: ${\Phi}$')
%% ------------------- Print the results -------------------------------
if nargout == 0
% these obtained system parameters will be used in numerical simulations
% Note that these parameters are all in the configuration space.
fprintf('\n\n\n\n\n\n\n\n\n\n\')
disp("====================== RESULTS FOR AMM DISCRETIZATION =================================")
disp("For Research Use: output-only modal analysis tool box v0.0")
disp("Copyright by Hewenxuan Li, hewenxuan_li@uri.edu")

format("short")
disp(['Time elapsed (sec):', num2str(time_taken)])

disp('The stiffness matrix:')
K
disp('The mass matrix:')
M
disp('The forcing weight:')
fw

format shortEng
disp('The natural frequencies (rad/sec)')
omega_n

disp("=======================================================================================")
else
fprintf('\n')
disp("================= FIXED-FIXED MEMBRANE AMM DISCRETIZATION COMPLETE! =========================")
disp("For Research Use: output-only modal analysis tool box v0.0")
disp("Copyright by Hewenxuan Li, hewenxuan_li@uri.edu")

format("short")
disp(['Time elapsed (sec):', num2str(time_taken)])
end
