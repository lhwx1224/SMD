function [M, K, T, Lambda, Phi, Phi_tilde, fw, omega_n] = AMM_uniform_cant_beam(Nam, dx, m0, EI0, L0, EIx_str, mx_str, Fx_str, FastGen)
% AMM_uniform_beam() discretize the Euler-Bernoulli beam into mass,
% stiffness matrices, and weighting to the foring function in configuration
% space using the caconical form of the uniform beam solution.
% 
% syntax: AMM_uniform_beam(m0, EI0, L0, EIx_str, mx_str, FastGen)
% 
% input: Nam - integer, number of assumed modes (Nam < 20 for this release)
%        dx - float, spatial resolution used to generate the assumed mode
%        shapes.
%        m0 - mass constant *(which is going to be applied to a function
%             that discribe the beam's mass distribution) 
%        EI0 - flexural rigidity constant *(which is going to be applied to
%        a function that discribe the beam's stiffness distribution) 
%        L0 - length of the beam
%        EIx_str - string, the expression of the flexural rigidity as a
%                 function of space variable x (use mass as 'm', flexural
%                 rigidity as 'EI', and length of the beam as 'L') 
%        mx_str - string, the expression of the mass density as a function
%                 of x (use mass as 'm', flexural rigidity as 'EI', and
%                 length of the beam as 'L')  
%        FastGen - string/numeric, flag of whether or not fast generation
%                  of roots to the trenscendental equation is used, the
%                  default is true (a set of predetermined roots up to the
%                  order of 19 is included) 
%
% output: M - mass matrix of the dicretized system, Nam-by-Nam
%         K - stiffness matrix of the dicretized system, Nam-by-Nam
%         T - modal transformation between the configuration space and the
%             modal space, Nam-by-Nam
%         Lambda - eigenvalues corresponding to the modal matrix, Nam-by-Nam
%         Phi - resolved modal basis functions phi(x), 1/dx-by-Nam
%         fw - forcing weight vector applied to the forcing vector
%         omega_n - natural frequencies of the discretized modes
%           
% NOTE: use mass as 'm', flexural rigidity as 'EI', and length of the beam
% as 'L', see example!
%
% example: 
% EIx_str = EI*(L - 0*x)/L;
% mx_str = m*(L - 0*x)/L;
% Fx = dirac(x - L0);
% m0 = 1;
% EI0 = 1;
% L0 = 1;
% [M, K, fw, omega_n] = AMM_uniform_beam(m0, EI0, L0, EIx_str, mx_str, FastGen)
%
%
% Copyright by Hewenxuan Li, Nov 2021
% Output-only modal analysis toolbox v0.0

%% UNIFORM EULER BERNOULLI BEAM AND ITS EIGENFUNCTIONS
if nargin < 9
    FastGen = 'True';
end
warning('off')
tic
syms m EI L
% Consider a uniform beam with constant mass and elastic constant
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
xs = single(0:dx:L0); % Sampled spatial vector

syms betaL x 
% Here, a change of variable is conducted to replace beta*L with betaL for
% solving the transcendental equation

objective = betaL*(1 + cos(betaL)*cosh(betaL))/(cos(betaL) + cosh(betaL)); % Equivalent to Eq, (o) pg. 401 of [1]
objs = double(subs(objective,betaL,[0:0.01:1000]));
if isequal(lower(FastGen),'true')
    % The center of the roots from solving the canonical form of the
    % objective function listed above
    loc_roots = [1.87510406871196,4.69409113297418,7.85475743823761,...
        10.9955407348755,14.1371683910465,17.2787595320882,20.4203522510413,...
    	23.5619449018064,26.7035375555183,29.8451302091028,32.9867228626928,...
    	36.1283155162826,39.2699081698724,42.4115008234622,45.5530934770520,...
    	48.6946861306418,51.8362787842316,54.9778714378214,58.1194640914112]';
else

% Peak finder 
beta_max = 1000;
disp('Finding roots of the designated problem')
disp(['Maximum beta value:', num2str(beta_max)])
xs_roots = 0:0.01:beta_max;
abs_d_obj = double(subs(abs(diff(objective,1)),betaL,xs_roots));
[pks, locs] = findpeaks(abs_d_obj);
loc_roots = xs_roots(locs)';
% ------------------------------------------------------------------------
plotrange = [0 beta_max];
clf
fplot(objective, [0 beta_max])
hold on
stem(loc_roots, max(pks)*ones(size(pks)), 'k--')
plot([plotrange(1) plotrange(end)], [0 0], 'k--')
xlabel('$\beta L$')
ylabel('Objective - $Y(\beta L)$')
legend('$Y(\beta L)$', '$\{\beta L|Y(\beta L) \approx 0\}$','location','southwest')
% ------------------------------------------------------------------------
end
disp(['Total number of modes solved:', num2str(length(pks))])
% This step solves for the variable z = beta*L from the objective
% transcendental equation. First, observe from the above plot and
% apporximately locate where are the roots located; then, proceed to
% numerical root finding. Otherwise, one can use the root-finding algorithm
% instead.

beta_range = [-0.5 0.5]; % threshold range for the root finder aside from the zeros
beta_hat = [loc_roots+beta_range(1), loc_roots+beta_range(2)];

% Define root vector
root = zeros(1,size(beta_hat,1));
% Find roots
for i = 1:length(root)
    root(i) = vpasolve(objective, betaL, beta_hat(i,:)); 
end

% Change the variable back to beta
beta = root/L0;
% --------------------- TESTING (to increase stability and accuracy) ------------
beta = sym(beta);
% ------------------------------------------
% Define the eigenfunction to the *UNIFORM BEAM* problem Y = Y(x) by
% plugging in the values solved from the transcendental equation.
Y = sin(beta*x) - sinh(beta*x) - ((sin(beta*L0) + sinh(beta*L0))./(cos(beta*L0) + cosh(beta*L0))).*(cos(beta*x) - cosh(beta*x));
% ddY = diff(Y,2);

% natural angular frequency (defined as numerical values)
% omega = (beta*L0).^2*sqrt(EI0/m0/L0^4); % See Example 8.4 Eq. (a)

% Calculate the magnitude of the eigenfunctions (via numerical int of the inner product of the eigenfunction)
Y_norm = sqrt(double(int(Y.^2,x,[0 L0])));

% Calculate the normalized eigenfunctions
Yn = Y./Y_norm;

% ---------------------- TESTING ----------------------
% Yns = zeros(length(xs),length(Yn));
% for i = 1:length(Y)
%     Yns(:,i) = double(subs(Yn(i), x, xs));
% end
% --------------------------------------------------
%% -------------------- Eigenvalue problem ----------------------
% syms EIx mx
EIx = str2sym(EIx_str);
mx = str2sym(mx_str);
Fx = str2sym(Fx_str);
EIx = subs(EIx, [EI L], [EI0 L0]);
mx = subs(mx, [m, L], [m0, L0]);
Fx = subs(Fx, L, L0);
% ----------------------- FIGURE 2 -------------------------------- 
EIs = double(subs(EIx, x, xs));
EIs = EIs.*ones(1,length(EIs));
ms = double(subs(mx, x, xs));
ms = ms.*ones(1,length(ms));
zline = zeros(1,length(EIs));
figure(2),clf
subplot(211)
vfill2curves(xs, EIs, zline), title('Flexural Rigidity Distribution'), ylabel('$EI(x)$')
ylim([0 EI0])
subplot(212)
vfill2curves(xs, ms, zline),title('Mass Distribution'),ylabel('$m(x)$'),xlabel('$x$')
sgtitle('Imposed $EI(x)$ and $m(x)$')
ylim([0 m0])
pause(0.01)
% Obtain the system matrices via numerical integration

% ----------------- Numerical Integration -------------------------------
% Initialize the integration for the stiffness matrix
K = zeros(Nam, Nam);
M = zeros(Nam, Nam);

for i = 1:Nam
    for j = 1:Nam
        fun_k = matlabFunction(EIx*diff(Yn(i), 2)*diff(Yn(j), 2));
        disp(['Calculating K(',num2str(i),',',num2str(j),')'])
        K(i, j) = integral(fun_k, 0, 1, "AbsTol",1e-3);
    end
end

round(K, 2);

% Initialize the integration for the mass matrix

for i = 1:Nam
    for j = 1:Nam
        fun_m = matlabFunction(mx*Yn(i)*Yn(j));
        disp(['Calculating M(',num2str(i),',',num2str(j),')'])
        M(i, j) = integral(fun_m, 0, 1, "AbsTol",1e-3);
    end
end

% -----------------------------------------------------------------------

% Solve for the system dynamics (Symmetric Eigenvalue Problem - EVP)
[evec, eval] = eig(inv(M)*K);
[eval_s, I] = sort(diag(eval));
eval = diag(eval_s);
evec = evec(:, I);
omega_n = sqrt(diag(eval));

% % -------------------- TEST CODE ----------------------------------
% % Obtain the modal mass and modal stiffness matrices after the assignments
% Km_am = evec'*K*evec;
% Mm_am = evec'*M*evec;
% mm_am = diag(Mm_am);                   % Extract modal mass matrix's diagonal elements
% km_am = diag(Km_am);                   % Extract modal stiffness matrix's diagonal elements
% 
% zeta_am = 0.01*(1:1:10)';
% 
% % Calculate the poles of the discretized system
% gamma_am = zeta_am.*omega_n;          % Real part of the pole
% omegad_am = omega_n.*sqrt(1 - zeta_am.^2);  % imaginary part of the pole
% 
% % Assign the modal damping matrix 
% Cm = diag(2*mm_am.*gamma_am);       % Modal damping matrix
% C = inv(evec')*Cm*inv(evec);
% A_am = [zeros(Nam), eye(Nam); -inv(M)*K, -inv(M)*C];
% [U_am, Lambda_am] = eig(A_am);
% [~, indx] = sort(abs(imag(diag(Lambda_am)))); % Sort the damped natural frequencies
% U_am = U_am(:,indx);
% % -----------------------------------------------------------------
% % Calculate the scaling to the new mode shapes
% Phi_x = Yn(1:Nam)*real(U_am(1:Nam, 1:2:end));
% Phi_x_norm = [];
% % Numerically integrate the norm of the new modal basis
%  for i = 1:Nam
%         fun_phi_x = matlabFunction(Phi_x(i).^2);
%         disp(['Calculating norm(Phi_',num2str(i),')'])
%         Phi_x_norm(i) = integral(fun_phi_x, 0, L0, "AbsTol",1e-3);
%  end
%  fplot(Phi_x,[0 1])
% % ---------------------- END TEST CODE ---------------------------------

% Sample the eigenfunction for a given resolution dx
Yns = zeros(length(xs), Nam);
for i = 1:size(Yns,2)
    Yns(:,i) = subs(Yn(i), x, xs);
end

% Impose the coordinate transformation provided by the EVP
% plot(Yns)

% Calculate the estimated mode shape for the non-uniform beam problem via
% modal transoformation
Phi = Yns*evec;
% ------------------- FIGURE 3 ------------------------------------
figure(3),clf
subplot(211)
plot(xs,Phi./Phi(end,:))
title('Non-uniform Beam via AMM (end-normalized)')
ylabel('Magnitude')
subplot(212)
plot(xs,Yns./Yns(end,:))
title('Uniform Beam Eigenfunctions')
ylabel('Magnitude')
xlabel('$x$')
%% Calculate the weights to the forcing function
% Define a force variable f(x,t) = F(x)G(t) and we only worries about the
% spatial part F(x) since the integration will be conducted over space.
syms f
fw = zeros(Nam, 1);
% Fx = dirac(x - L0); % <ANALYTICAL>
for i = 1:Nam
    disp(['Calculating fw(',num2str(i),')'])
    fw(i) = double(int(Fx*Yn(i),x,[0 L0])); % <ANALYTICAL>
end
% ------------------------------------------------------------------------
time_taken = toc;
%% Wrapping up the output variables
T = evec;
Lambda = eval;
Phi_tilde = Yns;
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
end
% ------------------------------------------------------------------------
