% Assumed-modes method (assuming the mode is at hand)
% References:
% [1] Fundamentals of Vibrations, Leonard Meirovitch, 1ed
% ------------------------------------------------------------------------
%
% This code is dedicated to the output-only modal analysis toolbox
% In this program, a uniform Euler-Bernoulli beam is considered as a
% reference to generate a set of eigenfunctions that are going to be used
% in assumed-mode series discretization (i.e., the Ritz basis functions).
% (1) a uniform EB beam is assumed and its eigenfunctions are resolved via
% numerical root finding. 
% (2) Then, the eigenfunctions are normalized and ready to be used as basis
% functions for assumed-modes modal analysis
% (3) The assumed modes (i.e., the eigenfunctions from the uniform beam)
% are used as basis to come up with the discretized mass, stiffness
% matrices, and forcing weight (namely, M, K, fw), respectively. Then, the
% mass, flexural rigidity, and forcing, as functions of space can be
% designated and used in the discretization. 
%
% Hewenxuan Li, Nov. 2021
%
% ------------------------------------------------------------------------

%% UNIFORM EULER BERNOULLI BEAM AND ITS EIGENFUNCTIONS
% Consider a uniform beam with constant mass and elastic constant
clear
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
m = 1;
EI = 1;
M = 1;
L = 1;    % Lets assume the length is one for now

syms betaL x 
% Here, a change of variable is conducted to replace beta*L with betaL for
% solving the transcendental equation

objective = betaL*(1 + cos(betaL)*cosh(betaL))/(cos(betaL) + cosh(betaL)); % Equivalent to Eq, (o) pg. 401 of [1]
% objective = cos(betaL)*cosh(betaL) + 1; % Eq, (o) pg. 401 of [1]
% f = (beta*(m + cosh(L*beta))*(m*cos(L*beta) - M*beta*sin(L*beta)) + M*beta*cos(L*beta)*sinh(L*beta))/(cos(L*beta) + cosh(L*beta));

% Peak finder 
xs_roots = 0:0.01:60;
abs_d_obj = double(subs(abs(diff(objective,1)),betaL,xs_roots));
[pks, locs] = findpeaks(abs_d_obj);
loc_roots = xs_roots(locs)';
% ------------------------------------------------------------------------
plotrange = [0 60];
clf
fplot(objective, [0 60])
hold on
stem(loc_roots, max(pks)*ones(size(pks)), 'k--')
plot([plotrange(1) plotrange(end)], [0 0], 'k--')
% plot(xs_roots, abs_d_obj)
% fplot(diff(objective,1), [0 60])
xlabel('$\beta L$')
ylabel('Objective - $Y(\beta L)$')
legend('$Y(\beta L)$', '$\{\beta L|Y(\beta L) \approx 0\}$','location','southwest')
% ------------------------------------------------------------------------

% This step solves for the variable z = beta*L from the objective
% transcendental equation. First, observe from the above plot and
% apporximately locate where are the roots located; then, proceed to
% numerical root finding. Otherwise, one can use the root-finding algorithm
% instead.

beta_range = [-0.5 0.5]; % threshold range for the root finder aside from the zeros
beta_hat = [loc_roots+beta_range(1), loc_roots+beta_range(2)];

% Estimated range of roots (if not default, end users have to modify)
% beta_hat = [1.8, 1.9; 4.5, 5; 7.6, 8.2; 10.8, 11.3; 13.9, 14.3; 17, 17.6; ...
%     20.2 20.5; 23.2, 23.8; 26.5, 26.9; 29.5, 29.9; 32.63, 33.34;...
%     35.9, 36.5; 39.1, 39.5; 42.2, 42.52; 45.4, 45.7; 48.6, 48.9; 51.6, 52;...
%     54.7, 55; 58, 58.3];

% Define root vector
root = zeros(1,size(beta_hat,1));
% Find roots
for i = 1:length(root)
    root(i) = vpasolve(objective, betaL, beta_hat(i,:)); 
end

% Change the variable back to beta
beta = root/L;

% Define the eigenfunction to the *UNIFORM BEAM* problem Y = Y(x)
Y = sin(beta*x) - sinh(beta*x) - ((sin(beta*L) + sinh(beta*L))./(cos(beta*L) + cosh(beta*L))).*(cos(beta*x) - cosh(beta*x));
ddY = diff(Y,2);

% Y_fun = matlabFunction(Y);
% ddY_fun = matlabFunction(ddY);

% ===== THIS PART OF CODE SHOWS THE EQUIVALENCE OF INT AND INTEGRAL =====
% % Compare the integration value between analytical integration and
% % numerical one
% qsym = double(int(Y(j),x,[0 1]));
% 
% fun = @(x, beta, L) sin(beta*x) - sinh(beta*x) - ((sin(beta*L) + sinh(beta*L))./(cos(beta*L) + cosh(beta*L))).*(cos(beta*x) - cosh(beta*x));
% qnum = integral(@(x) fun(x, beta(1), 1),0,1);
% qsym == qnum
% % -------------------------------------------------------------------------

% natural angular frequency (defined as numerical values)
omega = (beta*L).^2*sqrt(EI/m/L^4); % See Example 8.4 Eq. (a)

% Calculate the magnitude of the eigenfunctions (via numerical int of the inner product of the eigenfunction)
Y_norm = sqrt(double(int(Y.^2,x,[0 1])));

% Calculate the normalized eigenfunctions
Yn = Y./Y_norm;

% -------------------------- CODE FOR DEBUGGING ---------------------------
% Plot the obtained eigenfunction (default plot 10)
% figure(2)
% fplot(Y(:,1:10), [0 1])
% Plot the normalized eigenfunctions
% fplot(Yn(:,1:10), [0 1])
% -------------------------------------------------------------------------

%% -------------------- Eigenvalue problem ----------------------
syms EIx mx
EI0 = 1;
m0 = 1;
EIx = EI0*(L - 0*x)/L;
mx = m0*(L - 0.5*x)/L;

% ----------------------- FIGURE 2 -------------------------------- 
figure(2),clf
subplot(211)
fplot(EIx, [0 L]), title('Flexural Rigidity Distribution'), ylabel('$EI(x)$')
ylim([0 EI0])
subplot(212)
fplot(mx, [0 L]), title('Mass Distribution'),ylabel('$m(x)$'),xlabel('$x$')
sgtitle('Imposed $EI(x)$ and $m(x)$')
ylim([0 m0])
% EIx = EI0;
% mx = m0;
pause(0.01)
% Obtain the system matrices via numerical integration
Nam = 10;
% ---------------- Symbolic Integration ----------------------------------
% This part wont generally work for more complex x-dependent cases
% ------------------------------------------------------------------------
% Initialize the mass and stiffness matrices 
% Ka = zeros(Nam, Nam);
% Ma = zeros(Nam, Nam);
% % Initialize the integration for the sitffness matrix
% tic
% for i = 1:Nam
%     for j = 1:Nam
%         disp(['Calculating K(',num2str(i),',',num2str(j),')'])
%         Ka(i, j) = double(int(EIx*diff(Yn(i), 2)*diff(Yn(j), 2),x,[0 1]));
%     end
% end
% toc
% 
% round(Ka, 2);
% % Initialize the integration for the mass matrix
% tic
% for i = 1:Nam
%     for j = 1:Nam
%         disp(['Calculating M(',num2str(i),',',num2str(j),')'])
%         M(i, j) = double(int(mx*Yn(i)*Yn(j),x, 0, 1));
%     end
% end
% toc
% ------------------------------------------------------------------------

% ----------------- Numerical Integration -------------------------------
% % 
% Initialize the integration for the stiffness matrix
K = zeros(Nam, Nam);
M = zeros(Nam, Nam);

tic
for i = 1:Nam
    for j = 1:Nam
        fun_k = matlabFunction(EIx*diff(Yn(i), 2)*diff(Yn(j), 2));
        disp(['Calculating K(',num2str(i),',',num2str(j),')'])
        K(i, j) = integral(fun_k, 0, 1, "AbsTol",1e-3);
    end
end
toc
round(K, 2);

% Initialize the integration for the mass matrix
tic
for i = 1:Nam
    for j = 1:Nam
        fun_m = matlabFunction(mx*Yn(i)*Yn(j));
        disp(['Calculating M(',num2str(i),',',num2str(j),')'])
        M(i, j) = integral(fun_m, 0, 1, "AbsTol",1e-3);
    end
end
toc
% -----------------------------------------------------------------------

% Solve for the system dynamics (Symmetric Eigenvalue Problem - EVP)
[evec, eval] = eig(inv(M)*K);
[eval_s, I] = sort(diag(eval));
eval = diag(eval_s);
evec = evec(:, I);
omega_n = sqrt(diag(eval));
% Sample the eigenfunction for a given resolution dx
dx = 0.01;
xs = 0:dx:L;
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
Fx = dirac(x - L);

% Fx = 1;
for i = 1:Nam
    fw(i) = double(int(Fx*Yn(i),x,[0 L]));
end
% ------------------------------------------------------------------------
%% ------------------- Print the results -------------------------------
% these obtained system parameters will be used in numerical simulations
% Note that these parameters are all in the configuration space.
disp("====================== RESULTS FOR AMM DISCRETIZATION =================================")
disp("For Research Use: output-only modal analysis tool box v0.0")
disp("Copyright by Hewenxuan Li, hewenxuan_li@uri.edu")

format("short")
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

% ------------------------------------------------------------------------
