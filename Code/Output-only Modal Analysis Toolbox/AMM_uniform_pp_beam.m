function [M, K, T, Lambda, Phi, Phi_tilde, fw, omega_n] = AMM_uniform_pp_beam(Nam, dx, m0, EI0, L0, EIx_str, mx_str, Fx_str)
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
% Fx = dirac(x - L);
% m0 = 1;
% EI0 = 1;
% L0 = 1;
% [M, K, fw, omega_n] = AMM_uniform_beam(m0, EI0, L0, EIx_str, mx_str, FastGen)
%
%
% Copyright by Hewenxuan Li, Nov 2021
% Output-only modal analysis toolbox v0.0

%% UNIFORM EULER BERNOULLI BEAM AND ITS EIGENFUNCTIONS
disp("================= PIN-PINED BEAM AMM DISCRETIZATION STARTED! =========================")
addpath('L:\My Drive\Graduate study\Research\Projects\OS')
warning('off')
% Check if the parpool is running; if not, execute:
if isempty(gcp('nocreate'))
    clusterinfo = parcluster('local');
    Nworkers = clusterinfo.NumWorkers;
    parpool('local', Nworkers)
end
% MAIN starts here!
tic
syms m EI L betas x r 
% Consider a uniform beam with constant mass and elastic constant
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
xs = single(0:dx:L0); % Sampled spatial vector
 
% Here, a change of variable is conducted to replace beta*L with betaL for
% solving the transcendental equation
betas = r*pi/L0;

beta_r = subs(betas,r,1:1:Nam);

disp(['Total number of modes considered:', num2str(Nam)])
% This step solves for the variable z = beta*L from the objective
% transcendental equation. First, observe from the above plot and
% apporximately locate where are the roots located; then, proceed to
% numerical root finding. Otherwise, one can use the root-finding algorithm
% instead.

% Define the eigenfunction to the *UNIFORM PIN-PINED BEAM* problem Y = Y(x) by
% plugging in the values solved from the transcendental equation.
Y = sin(beta_r*x);

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
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])
% Obtain the system matrices via numerical integration

% ----------------- Numerical Integration -------------------------------
% Initialize the integration for the stiffness matrix
K = zeros(Nam, Nam);
M = zeros(Nam, Nam);

% count = 1;
% tic
% for i = 1:Nam
%     for j = 1:Nam
%         fun_k = matlabFunction(EIx*diff(Yn(i), 2)*diff(Yn(j), 2));
%         K(i, j) = integral(fun_k, 0, 1, "AbsTol",1e-3);
%         progress_bar((i-1)*Nam+j,Nam*Nam,'Caculating Stiffness Matrix')
%         count = count + 1;
%     end
% end
% toc
% Initialize the integration for the mass matrix
% count = 1;
% for i = 1:Nam
%     for j = 1:Nam
%         fun_m = matlabFunction(mx*Yn(i)*Yn(j));
%         M(i, j) = integral(fun_m, 0, 1, "AbsTol",1e-3);
%         progress_bar((i-1)*Nam+j,Nam*Nam,'Caculating Mass Matrix')
%         count = count + 1;
%     end
% end
% ----------------------------------------------------------------

% -------------------- TEST PARALLELIZATION ----------------------
% During the comparison of the single thread calc and the multicore calc
% (with 100 modes involved), the single core spent 322.37 secs compared to
% 80 secs obtained from the parallelization. This result suggests using
% parallelization whenever a large scale problem is involved (high modal
% structure simulation).
% Such results can be also obtained when only 10 modes involved: single
% core yielded 4.18 secs and multicore yielded 0.86 secs.

disp('Calculating Stiffness Matrix (Parallelized)')
parfor i = 1:Nam
    for j = 1:Nam
        fun_k = matlabFunction(EIx*diff(Yn(i), 2)*diff(Yn(j), 2));
        K(i, j) = integral(fun_k, 0, 1, "AbsTol",1e-3);
    end
end

disp('Calculating Mass Matrix (Parallelized)')
parfor i = 1:Nam
    for j = 1:Nam
        fun_m = matlabFunction(mx*Yn(i)*Yn(j));
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

% Sample the eigenfunction for a given resolution dx
% IF single core: 208.66 secs; if multicore: 4.15 secs
disp('Calculating Forcing Weights (Parallelized)')
Yns = zeros(length(xs), Nam);

parfor i = 1:Nam
    Yns(:,i) = subs(Yn(i), xs);

%     Yns(:,i) = subs(Yn(i), x, xs);
%     progress_bar(i,Nam,'Generating Sampled Mode Shapes')
end

% Impose the coordinate transformation provided by the EVP
% plot(Yns)

% Calculate the estimated mode shape for the non-uniform beam problem via
% modal transoformation
Phi = Yns*evec;

%% Calculate the weights to the forcing function
% Define a force variable f(x,t) = F(x)G(t) and we only worries about the
% spatial part F(x) since the integration will be conducted over space.
syms f
fw = zeros(Nam, 1);
% Fx = dirac(x - L0); % <ANALYTICAL>

for i = 1:Nam
    fw(i) = double(subs(int(Fx*Yn(i),x,[0 L0]), L, L0)); % <ANALYTICAL>
    progress_bar(i,Nam,'Calculating Forcing Weights')
end

% ------------------------------------------------------------------------
time_taken = toc;
%% Wrapping up the output variables
T = evec;
Lambda = eval;
Phi_tilde = Yns;

% ------------------- FIGURE 3 ------------------------------------
Indx10 = round(linspace(1,Nam,10));
figure(3),clf
subplot(411)
plot(xs,Phi(:, Indx10)./max(Phi(:, Indx10)))
title('Non-uniform Beam via AMM (end-normalized) ${\Phi}$')
ylabel('Magnitude')
subplot(412)
plot(xs,Yns(:, Indx10)./max(Yns(:, Indx10)))
title('Uniform Beam Eigenfunctions $\widetilde\Phi$')
ylabel('Magnitude')
xlabel('$x$')
subplot(413)
plot(fw)
xlabel('Mode $\#$')
ylabel('Forcing Weight')
subplot(414)
plot(abs(diag(evec)))
xlabel('Mode $\#$')
ylabel('$\mathrm{abs}(\mathrm{diag}(T))$')
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])
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
fprintf('\n\n\n\n\n\n\n\n\n\n\')
disp("================= PIN-PINED BEAM AMM DISCRETIZATION COMPLETE! =========================")
disp("For Research Use: output-only modal analysis tool box v0.0")
disp("Copyright by Hewenxuan Li, hewenxuan_li@uri.edu")

format("short")
disp(['Time elapsed (sec):', num2str(time_taken)])
end
addpath('L:\My Drive\Graduate study\Research\Projects\OS')
finishbeep('F')
% ------------------------------------------------------------------------
