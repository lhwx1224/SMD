% Forced Beam Assumed Modes Method (AMM) - FIXED-FIXED MEMBRANE SIMULATOR
% Simulation Environment Developed for the Modal Identification Toolbox
% Cantilever beam simulator using the assumed mode method and a wide variaty of forcing. This package can also save the data for later use. 
% Features:
% (1) Forcing type: Free Decay, Harmonic, Random (psudo-and-lowpass-filtered), Burst-random, Impulse;
% (2) Forcing Location: Beam and Base
% (3) Initial Excitation: various weight selection to ensure only particular modes involved in the response
% (4) Simulation Method: Numerical ODE solver v.s. Linear system iteration through transfer function/ polynomial filtering 
% (5) Sampling criteria and spatial resolution and more.
%  Copyright by Hewenxuan Li, April 2021.

% On-demand Reset
clear
close all
% Environment Setup
set(0,'DefaultFigureWindowStyle','docked')
DefaultPath = split(cd,'\');
for i = 1:length(DefaultPath)
    if isequal(DefaultPath{i},'SMD')
        flag = i;
        break
    end
end
DefaultPath = strjoin(DefaultPath(1:flag),'\');
addpath([DefaultPath,'\Code\OS'])
addpath([DefaultPath,'\Code\Output-only Modal Analysis Toolbox'])
set(0,'DefaultFigureWindowStyle','docked')
[Markers, LineStyles, Mycolors] = mystyle();   % Loading printing and plotting style
ReadSystem = 0;                                % Read System Matrix
%%
if ~ReadSystem

% ------------------------------------------------------------------------
%                      General System Setup
% ------------------------------------------------------------------------
Method = 'LSIM';         % Simulation Method
SimSpace = 'Modal';      % Simulation runs in which space
DampingType = 'Modal';   % Type of damping involved in the simulation

% ------------------------------------------------------------------------
%                Forcing Setup for the Simulation
% ------------------------------------------------------------------------
Forcing = 'FreeDecay';   % External Excitation Form
ExciteForm = 'Beam';     % Excitation Formulation (Beam - original, Base - d^2(original)/dt^2)
ForceWeight = 'Base';    % Weight used in weighting the force in the AMM
ExciteLoc = 10;          % Excitation location

% ------------------------------------------------------------------------
%                 Beam Mode Shapes Setup   
% ------------------------------------------------------------------------
% BEAM SPATIAL INFOMRATION
mm = 5;            % Number of eigenfunctions used in the x direction
nn = 2;            % Number of eigenfunctions used in the y direction
nm = mm*nn;        % Total number of modes considered in the simulation
a0 = 1;            % Length of the side in the x direction
b0 = 2;            % Length of the side in the y direction
% Membrane parameters (defuault mass density and the tension)
rho0 = 1;                            % Mass density constant
T0 = 1;                              % Tension function
% Mass density distribution function
rho_xy_str = '2*rho - (x - a/2).^2 - (y - b/2).^2';
% Flexural rigidity distribution
T_xy_str = 'T*heaviside(x)*heaviside(y)';                
% Designate the forcing function distribution
if isequal(lower(ForceWeight), 'beam')
    F_xy_str = strcat('dirac(x - L*',num2str(ExciteLoc),'/10)');
elseif isequal(lower(ForceWeight), 'base')
    F_xy_str = '1';
end
% 
dx = 0.01*a0;         % Resolution in the x direction
dy = 0.01*b0;         % Resolution in the y direction

x = 0:dx:a0;                          % Define domain of the membrane
y = 0:dy:b0;
p = length(x);                        % # of sampling points from both ends
% SIMULATION TIME SETUP
dt = 0.05;                        % sampling time (Default for nm == 10 is 0.0025)
fs = 1/dt;                          % sampling frequency
tend = 60;                           % ending time of simulation
t = 0:dt:tend;                      % Time vector
rr = 1;                             % set the resapling rate (optional). Test other than random should use rr = 1;

% EIGENVALUE PROBLEM VIA AMM_uniform_beam FUNCTION
[M_am, K_am, T, Lambda, phi_am, phi_tilde, weight, omega_am] = AMM_uniform_FF_membrane(mm, nn, dx, dy, a0, b0, rho0, T0, rho_xy_str, T_xy_str, F_xy_str);

% Plot basic information of the simulated system
plot(omega_am*2*pi,'.')
pbaspect([1.618,1,1])
grid on
xlabel('Mode Index')
ylabel('Natural Frequency (Hz)')

% SAVE THE SYSTEM INFORMATION FOR FUTURE USE
run SaveSystemInfo_FFM.m

else
   load('System.mat');
    M_am = System.M;
    K_am = System.K;
    T = System.T;
    Lambda = System.Lambda;
    phi_am = System.phi_am;
    phi_tilde = System.phi_tilde;
    weight = System.weight;
    omega_am = System.omegan;
    tend = System.tend;
    fs = System.fs;
    dt = System.dt;
    dx = System.SpatialResolution.dx;
    dy = System.SpatialResolution.dy;
    a0 = System.Dimensions.a0;
    b0 = System.Dimensions.b0;
    rho0 = System.UniformDensity;
    T0 = System.UniformTension;
    rr = System.rr;
    % Number of modes for the discretized model
    nm = System.nm;                             % Number of modes to consider
    mm = System.mm;
    nn = System.nn;
    rho_xy_str = System.Density;
    T_xy_str = System.Tension;
    F_xy_str = System.Force;
    Forcing = System.ForcingType;
    Method = System.SimulationMethod;
    SimSpace = System.SimulationSpace;
    ExciteForm = System.ExciteForm;
    ForceWeight = System.ForceWeight;
    DampingType = System.DampingType;
    ExciteLoc = System.ExciteLocation; % Excitation location
    t = 0:dt:tend;                      % Time vector
    x = 0:dx:a0;                          % Define domain of the beam
    y = 0:dy:b0;
    p = length(x);
end
%% EXTRACT COMPONENTS OF THE NATURAL FREQUENCIES

% ------------------------------------------------------------------------
%            MODAL ANALYSIS SECTION STARTS
% NOTE: modal analysis only in terms of solving the modal dynamics; the
% results will be later used as function $q(t) = T^{T}*\eta(t)$.
% ------------------------------------------------------------------------

% Recast the mode tensor (in cell format) into matrices
% Phi_tilde = zeros((length(x)-2)*(length(y)-2), mm*nn);  % allocate memory for assumed mode shapes
% Phi_am = zeros((length(x)-2)*(length(y)-2), mm*nn);  % allocate memory for approximated mode 
Phi_tilde = zeros((length(x))*(length(y)), mm*nn);  % allocate memory for assumed mode shapes
Phi_am = zeros((length(x))*(length(y)), mm*nn);  % allocate memory for approximated mode shapes
for i = 1:length(phi_tilde)
    % Remove the boundaries for both Phi_am and Phi_tilde
    phi_tilde_temp = phi_tilde{i}; % place holder for phi_tilde{i}
%     phi_tilde_temp(:,[1 size(phi_tilde_temp,2)]) = []; % remove the first and last columns 
%     phi_tilde_temp([1 size(phi_tilde_temp,1)],:) = []; % remove the first and last rows
    Phi_tilde(:,i) = reshape(phi_tilde_temp, [], 1);   % reshape it into a column vector and store in matrix Phi_tilde
    % Remove the boundaries for both Phi_am and Phi_am
    phi_am_temp = phi_am{i}; % place holder for phi_tilde{i}
%     phi_am_temp(:,[1 size(phi_am_temp,2)]) = []; % remove the first and last columns 
%     phi_am_temp([1 size(phi_am_temp,1)],:) = []; % remove the first and last rows
    Phi_am(:,i) = reshape(phi_am_temp, [], 1);   % reshape it into a column vector and store in matrix Phi_tilde
end

% Build Matrices (Dynamic matrices in configuration space)
n = nm;                          % number of total modes

% Obtain the modal mass and modal stiffness matrices after the assignments
Km_am = T'*K_am*T;
Mm_am = T'*M_am*T;
mm_am = diag(Mm_am);             % Extract modal mass matrix's diagonal elements
km_am = diag(Km_am);             % Extract modal stiffness matrix's diagonal elements
Km_am = diag(km_am);             % Only retain the diagonal terms
Mm_am = diag(mm_am);

%% ASSIGN DAMPING MATRIX
close all
% ------------- CONDUCT THE ASSIGNMENT IN THE MODAL SPACE ----------------
% Assign modal damping coefficients
zeta_am0 = 0.01;
zeta_am = zeta_am0*ones(nm,1);

% Calculate the poles of the discretized system
gamma_am = zeta_am.*omega_am;          % Real part of the pole
omegad_am = omega_am.*sqrt(1 - zeta_am.^2);  % imaginary part of the pole

% Assign the modal damping matrix 
Cm_am = diag(2*mm_am.*gamma_am);       % Modal damping matrix

if isequal(lower(DampingType), 'rayleigh')
% Determine the Rayleigh damping coefficients (according to the ith and jth modal damping ratios)
alpha = 2*omega_am(1)*omega_am(2)/(omega_am(2)^2 - omega_am(1)^2)*[omega_am(2), -omega_am(1); -1/omega_am(2), 1/omega_am(1)]*[zeta_am(1); zeta_am(2)];
C_am = alpha(1)*M_am + alpha(2)*K_am;   % Assign Rayleigh Damping
cm_am = diag(T'*C_am*T);                % Modal damping coefficients
zeta_am = cm_am./(2*mm_am.*omega_am);   % Re-calculate the damping ratios
elseif isequal(lower(DampingType), 'modal')
    C_am = inv(T')*Cm_am*inv(T);
end
% Add a damper in the middle of the beam
C0 = 0; % damping coefficient of the damper in the middle
iC0 = 5; % location of the lumped damper
% Append non-proportional damping effect to the damping matrix
C_am = C_am + C0*Phi_am(iC0, :)'*Phi_am(iC0,:);

% ------------------------------------------------------------------------
% State-space Dynamic Matrix in Configuration Space*
% ------------------------------------------------------------------------
% A_am = [zeros(n), Mm_am; -Km_am, -Cm_am]; 
A_am = [zeros(n), eye(n); -inv(M_am)*K_am, -inv(M_am)*C_am];
B_am = [zeros(n); inv(M_am)];
[U_am, Lambda_am] = eig(A_am);
[U_am_T, ~] = eig(A_am');

% For the right eigenvectors
[omegads, indx] = sort(abs(imag(diag(Lambda_am))), 'ascend'); % Sort the damped natural frequencies
omegads = omegads(1:2:end);   % extract the damped natural frequencies of the modes
fds = omegads/2/pi;
zetas = - real(diag(Lambda_am)); % poles (complex eigenvalues)
zetas = zetas(indx);          % sort the real part accroding to the damped natural frequencies
zetas = zetas(1:2:end);       % Extract the unique damped natural frequencies
omegaeig = abs(diag(Lambda_am));  % Extract the undamped natural frequencies
omegaeig = omegaeig(indx);    % Sort the extracted undamped natural frequencies
omegaeig = omegaeig(1:2:end);
zetas = zetas./omegaeig;      % Extract the modal dampling ratios 
cmodes = [Phi_tilde*U_am(1:n, indx); Phi_tilde*U_am(n+1:2*n, indx)]; % Modes for the actual physical space <-> modal space
rmodes = real(cmodes(1:size(Phi_tilde,1),1:2:end)); % Real part of the modes 
imodes = imag(cmodes(1:size(Phi_tilde,1),1:2:end)); % Imaginary part of the modes
% rmodesn = [zeros(1,size(rmodes,2)); rmodes/diag(max(rmodes))]; % zero-padding and normalization
% imodesn = [zeros(1,size(imodes,2)); imodes/diag(max(imodes))];

% Put the mode shapes into order
U_am = U_am(:,indx);
U_am_T = U_am_T(:,indx);

U_temp = zeros(size(U_am));
U_T_temp = zeros(size(U_am_T));
for i = 1:size(U_temp,2)
    if mod(i, 2) == 0
        U_temp(:,i) = U_am(:,i-1);
        U_T_temp(:,i) = U_am_T(:,i-1);
    else
        U_temp(:,i) = U_am(:,i+1);
        U_T_temp(:,i) = U_am_T(:,i+1);
    end
end

lambda_am = diag(Lambda_am);
Lambda_am = diag(lambda_am(indx));

% For the adjoint modes (Left eigenvectors)
% cmodesa = [Phi_tilde*U_am_T(1:n, indx); Phi_tilde*U_am_T(n+1:2*n, indx)]; % Modes for the actual physical space <-> modal space
% rmodesa = real(cmodes(1:size(Phi_tilde,1),1:2:end)); % Real part of the modes 
% imodesa = imag(cmodes(1:size(Phi_tilde,1),1:2:end)); % Imaginary part of the modes
% rmodesa = [zeros(1,size(rmodesa,2)); rmodesa/diag(rmodesa(end,:))]; % zero-padding and normalization
% imodesa = [zeros(1,size(imodesa,2)); imodesa/diag(imodesa(end,:))];
% 
Phi_x = (Phi_tilde*U_am(1:n, 1:2:end)); 
Phi_v = (Phi_tilde*U_am(n+1:end, 1:2:end));
for i = 1:n
    norm_Phi_x(i) = norm(Phi_x(:,i));
end

% ------------------------------------------------------------------------
% State-space dynamic matrix in the modal space (psudo-modal space)
% ------------------------------------------------------------------------
figure
plot(fds)
xlabel('Mode $\#$')
ylabel('$f_d$')
grid on
set(gcf,'papersize',[6 2.5]),set(gcf,'paperposition',[0 0 6 2.5]);
% Define the modal state-space system matrix
Am_am = [zeros(n), eye(n); -inv(Mm_am)*Km_am, -inv(Mm_am)*Cm_am];
% Define the modal controllability matrix
Bm_am = [zeros(n); inv(Mm_am)*T'];
% Check if the system is proportionally damped:
if sum(round(K_am*inv(M_am)*C_am, 4) == round(C_am*inv(M_am)*K_am,4),"all") == n^2
    disp('------ System is PROPORTIONALLY damped! ------')
    npd = 0;
else
    disp('====== System is NON-PROPORTIONALLY damped! ======')
    npd = 1;
end
%% ===================== Get forcing function =============================
switch lower(Forcing)
    % -------------------------- HARMONIC FORCING -----------------------
    case 'harmonic'
        loading = "Harmonic";                              % Set the loading type
        rr = 1;                                            % Force the rr to be one 
        timing.Tstart = 0;                                 % Starting time of simulations
        timing.Tend = tend;                                % Ending time of simulations
        timing.dT = dt*rr;                                 % Sampling time interval
        finfo.A = 2;                    % Input amplitude
        finfo.omega = 15;                % Input frequency
        finfo.windtype = [];                               % Window type
        finfo.pwind = [];                                  % Window size
        finfo.floc = lower(ExciteForm);
        Tend = timing.Tend;
        % Generate forcing
        [forcing, ~, tspan, fs, dt, ~] = forcingfnc(loading, timing, finfo, rr); 
        disp('INSPECTION OF THE FORCING FUNCTION:')
        DataInspection(ppval(forcing,tspan),tspan,fs)
        pause(0.1)
    % -------------------------- IMPULSE FORCING -----------------------
    case 'impulse'
        loading = "Impulse";            % Set the loading type
        rr = 1;                         % Force the rr to be one 
        timing.Tstart = 0;              % Starting time of simulations
        timing.Tend = tend;             % Ending time of simulations
        timing.dT = dt/rr;              % Sampling time interval
        finfo.A = 2;                    % Input amplitude (area under pulse)
        finfo.omega = [];               % Input frequency
        finfo.windtype = "Harmonic";    % Window type
        finfo.pwind = 2;                % Window size (pulse width = pwind*dt)
        finfo.Nimp = 1;
        finfo.floc = lower(ExciteForm);
        Tend = timing.Tend;
        [forcing, ~, tspan, fs, dt, ~] = forcingfnc(loading, timing, finfo, rr);
        disp('INSPECTION OF THE FORCING FUNCTION:')
        DataInspection(ppval(forcing,tspan),tspan,fs)
        pause(0.1)
    % -------------------------- RANDOM FORCING -------------------------
    case 'random'
        loading = "Random";
        timing.Tstart = 0;
        timing.Tend = tend;
        timing.dT = dt*rr;
        finfo.A = 10;
        finfo.omega = [];
        finfo.windtype = "Rectangular";
        finfo.pwind = 0;
        finfo.floc = lower(ExciteForm);
        Tend = timing.Tend;
        [forcing, ~, tspan, fs, dt, ~] = forcingfnc(loading, timing, finfo, rr);
        disp('INSPECTION OF THE FORCING FUNCTION:')
        DataInspection(ppval(forcing,tspan),tspan,fs)
        pause(0.1)
     % -------------------------- BURST-RANDOM FORCING ------------------
     case 'burstrandom'
        loading = "Random";
        timing.Tstart = 0;
        timing.Tend = tend;
        timing.dT = dt*rr;
        finfo.A = 1;
        finfo.omega = [];
        finfo.windtype = "BurstRandom";
        finfo.N = 4;
        finfo.pwind = [];
        finfo.floc = lower(ExciteForm);
        Tend = timing.Tend;
        [forcing, ~, tspan, fs, dt, ~] = forcingfnc(loading, timing, finfo, rr);
        disp('INSPECTION OF THE FORCING FUNCTION:')
        DataInspection(ppval(forcing,tspan),tspan,fs)
        pause(0.1)
    % -------------------------- FREE RESPONSE --------------------------
    case 'freedecay'
        loading = "FreeDecay";
        timing.Tstart = 0;
        timing.Tend = tend;
        timing.dT = dt*rr;
        finfo.A = 0;
        finfo.omega = [];
        finfo.windtype = [];
        finfo.N = [];
        finfo.pwind = [];
        finfo.floc = lower(ExciteForm);
        Tend = timing.Tend;
        [forcing, ~, tspan, fs, dt, ~] = forcingfnc(loading, timing, finfo, rr);
        disp('INSPECTION OF THE FORCING FUNCTION:')
        DataInspection(ppval(forcing,tspan),tspan,fs ...
            )
        pause(0.1)
end

% Display the adjusted sampling frequency (with resampling ratios ~= 1)
disp(['ADJUSTED SAMPLING FREQUENCY:', num2str(fs),'Hz'])
%% RUN SIMULATIONS
% ------------------------------------------------------------------------
%                  SIMULATION STARTS HERE
% ------------------------------------------------------------------------
close all
% Initial value problem setup
IW = zeros(10,1);                      % Initial weight vector
InitialWeight = 1;                     % Initial weight value
q0indx = 1:10;         % To which mode the weight is assigned
IW(q0indx) = InitialWeight;
% ASSIGN THE MODAL PARTICIPATION FOR INITIAL VALUE PROBLEM
IW = linspace(1000,1,nm);
IW = exp(-0.1*(1:1:nm));
% IW = ones(1,nm);
IWindx1 = 1:2:2*nm;
IWindx2 = 2:2:2*nm;
IW_temp = [[IWindx1; IW],[IWindx2; IW]];
[~, I] = sort(IW_temp(1,:));
IW = IW_temp(:,I);
% IW(1:2) = [0.01 0.01];

% MODES INCLUDED
ModePool = [1:2:1000; 2:2:1000];                    % Generate a large enough mode index pool
ModePool = [1:1:size(ModePool,2); ModePool];
ModeIndices = [1:nm];                                  % Input the included mode indices
IncModes = reshape(ModePool(2:3,ModeIndices),[],1); % The included complex mode shapes
IncModes = sort(IncModes, 'ascend');                % Sort the complex modal index matrix
% IMPOSE SAME INITIAL CONDITION WHICH IS INDPT OF MODES
xs = round(linspace(1, p - 1, nm));
% e0 = [Phi_am(xs,1)];  % INITIAL DISPLACEMENT CONDITION
TUx_am = T*U_am(1:n,:);
q0 = real(sum(U_am(1:n,IncModes)*diag(IW(2,IncModes)),2));  % INITIAL DISPLACEMENT

% Simulate the dynamical system in the modal space
if isequal(lower(SimSpace), 'modal')
    A = Lambda_am;          % IF modal simulation, A is the modal dynamic matrix
    B = inv(U_am)*B_am;     % B is the modally transformed B_am
elseif isequal(lower(SimSpace), 'config')
    A = A_am;               % IF configuration simulation, A = A_am; B = B_am
    B = B_am;
elseif isequal(lower(SimSpace), 'pmodal')
    A = Am_am;              % IF Psudo-modal simulation, A = configuration space modal matrix
    B = Bm_am;              
end

if isequal(lower(Forcing),'freedecay')
    if isequal(lower(SimSpace), 'pmodal')
    eta0 = [q0; zeros(n,1)];            % Initial conditions  
    elseif isequal(lower(SimSpace), 'modal')
    % Assign modal-space initial condition
    eta0 = inv(U_am)*[q0; zeros(n, 1)];
    elseif isequal(lower(SimSpace), 'config')
    % Assign configuration space initial condition
    eta0 = [q0; zeros(n,1)];            % Initial conditions  
    % ----------------------------------------------
    end
else
    eta0 = zeros(2*n,1);    % Initial conditions (default zeros(2*n,1))   
end

% ========================================================================
%       NUMERICAL ITERATION OR STEPPING FOR THE TIME COORDINATES
% ========================================================================
switch lower(Method)
    case 'lsim'
 % ------------------- Linear System Simulation "lsim" -------------------
        Fm = repmat(ppval(forcing,tspan),nm,1)'; % Augmented forcing matrix
        Um = zeros(size(Fm));                    % allocate weighted forcing matrix
        for i = 1:length(weight) % The for loop can be replaced by a repmat of weight
            Um(:,i) = weight(i)*Fm(:,i); % weighted augmented forcing matrix
        end
        sys = ss(A, B, eye(size(A)), zeros(size(B))); % Define state space system in transfer function form
        tic
        [eta_am, t] = lsim(sys, Um, tspan, eta0);                      % Solve theplo system
        toc
        disp('INSPECTION OF THE MODAL COORDINATES:')
        if isequal(lower(SimSpace), 'pmodal')
            DataInspection(real(eta_am(:,1:nm)),t,fs)
        else
            DataInspection(real(eta_am),t,fs)
            sgtitle('Real Part of Modal Motion')
            DataInspection(imag(eta_am),t,fs)
            sgtitle('Imaginary Part of Modal Motion')
        end
        if isequal(lower(SimSpace), 'modal')
            y_am = eta_am(:,1:1:end)*U_am(1:n,1:1:end).'*Phi_tilde'; %<Modal - Config>
            dy_am = eta_am(:,1:1:end)*U_am(n+1:2*n,1:1:end).'*Phi_tilde'; %<Modal - Config>
            y_am = real(y_am);
            dy_am = real(dy_am);
        elseif isequal(lower(SimSpace), 'pmodal')
            y_am = eta_am(:,1:n)*T'*Phi_tilde';
            dy_am = eta_am(:,n+1:2*n)*T'*Phi_tilde';
        end

% -------------------- Simulation via Numerical Integration --------------
    case 'ode'
        % Numerical integration
        odefcn = @(t,q) GenLinBeam_conf(t,q,A,B,forcing,weight); % DEFINE ODE FUNCTION
        % options = odeset('OutputFcn',@(t,y,flag) odeoutput_bw(t,y,flag,m,a,k,B,alpha,beta1,beta2,gamma,c,n,f));
        disp('Numerical integration started:')
        tic
        [t,eta_am] = ode45(odefcn,tspan,eta0);
        toc
        disp('Numerical integration finished!')
        disp('INSPECTION OF THE MODAL COORDINATES:')
        if isequal(lower(SimSpace), 'modal')
            DataInspection(real(eta_am(:,1:nm)),t,fs)
        else
            DataInspection(real(eta_am),t,fs)
            sgtitle('Real Part of Modal Motion')
            DataInspection(imag(eta_am),t,fs)
            sgtitle('Imaginary Part of Modal Motion')
        end
        if isequal(lower(SimSpace), 'modal')
            y_am = eta_am(:,1:1:end)*U_am(1:n,1:1:end).'*Phi_tilde'; %<Modal - Config>
            dy_am = eta_am(:,1:1:end)*U_am(n+1:2*n,1:1:end).'*Phi_tilde'; %<Modal - Config>
            y_am = real(y_am);
            dy_am = real(dy_am);
        elseif isequal(lower(SimSpace), 'pmodal')
            y_am = eta_am(:,1:n)*T'*Phi_tilde';
            dy_am = eta_am(:,n+1:2*n)*T'*Phi_tilde';
        end
end

disp('INSPECTION OF THE PHYSICAL COORDINATES:')
DataInspection(real(y_am),t,fs)
sgtitle('Physical Space Displacement Response')
DataInspection(real(dy_am),t,fs)
sgtitle('Physical Space Velocity Response')

figure(5),clf
subplot(311)
yyaxis left
p1 = plot(real(Phi_x(:,ModeIndices(1:2:6))));
hold on
yyaxis right
p2 = plot((real(y_am(1,:))));
legend([p1(1) p2],{'$q_i(0)\phi_i,\,i = 1,2,3$','$\sum_{i=1}^{N}q_i(0)\phi_i$'})
axis tight
pbaspect([2.5 1 1])
subplot(312)
yyaxis left
plot(abs(eta_am(1,1:2:size(eta_am,2))),'o')
set(gca,'yscale','log')
ylabel('Initial Modal Amplitude --- $q_i(0)$')
yyaxis right
plot(var(eta_am(:,1:2:size(eta_am,2))),'o')
ylabel('Var. of $q_i(t)$ --- $E[q_i(t)]^2$')
axis tight
xlabel('Mode Index')
pbaspect([2.5 1 1])
subplot(313)
imagesc(real(y_am'))
xlabel('Time Sample')
ylabel('Spatial Sample')
mycolorbar()
pbaspect([2.5 1 1])
set(gcf,'papersize',[6 6*1.5])
set(gcf,'paperposition', [0 0 6  6*1.5])
print('ICnField_ff_membrane_1.png','-dpng','-r600')

y_modal = y_am;      % Response matrix
dy_modal = dy_am;    % Velocity response matrix

% Convert 2-d array of system response into a 3-d tensor
% y_modal_tensor = zeros(length(x) - 2, length(y) - 2, size(y_modal,1));
y_modal_tensor = zeros(length(x), length(y), size(y_modal,1));
dely_modal_tensor = zeros(length(x), length(y), size(y_modal,1));
dely_modal = zeros(size(y_modal));
for i = 1:size(y_am,1)
%     y_modal_tensor(:,:,i) = reshape(y_modal(i,:), length(x) - 2, length(y) - 2);
    y_modal_tensor(:,:,i) = reshape(y_modal(i,:), length(x), length(y));
    [dely_modal_tensor_temp, ~] = GenFiniteDiff2(y_modal_tensor(:,:,i), dx, dy);
    dely_modal_tensor(:,:,i) = dely_modal_tensor_temp;
    dely_modal(i,:) = reshape(dely_modal_tensor(:,:,i), length(x)*length(y),1);
end
%% VISUALIZING CONFIGURATION SPACE RESPONSE
close all

% Find the mininum and maximum values in the tensor
min_y_tensor = min(y_modal_tensor, [], 'all');
max_y_tensor = max(y_modal_tensor, [], 'all');
min_dy_tensor = min(dely_modal_tensor, [], 'all');
max_dy_tensor = max(dely_modal_tensor, [], 'all');

% Ticklabels
xtickvec = 0:round(length(x)/5):length(x);
ytickvec = 0:round(length(y)/5):length(y);
xvec = linspace(0, a0, length(xtickvec));
yvec = linspace(0, b0, length(ytickvec));

figure(6),clf
for i = 1:2:size(y_am,1)
    subplot(121)
    s = surf(y_modal_tensor(:,:,i));
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    xticks(xtickvec);
    yticks(ytickvec);
    xticklabels(cellstr(string(xvec)));
    yticklabels(cellstr(string(yvec)));
    zlim([min_y_tensor max_y_tensor])
    xlim([1 length(x)])
    ylim([1 length(y)])
    mycolorbar();
    caxis([min_y_tensor max_y_tensor])
    view(-30, 30)
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$w(x,y)$')
    title(['Frame:', num2str(i), ' out of ',num2str(size(y_am,1))])
    pbaspect([a0 b0 0.5*a0])
    subplot(122)
    s = surf(dely_modal_tensor(:,:,i));
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    xticks(xtickvec);
    yticks(ytickvec);
    xticklabels(cellstr(string(xvec)));
    yticklabels(cellstr(string(yvec)));
    zlim([min_dy_tensor max_dy_tensor])
    xlim([1 length(x)])
    ylim([1 length(y)])
    mycolorbar();
    caxis([min_dy_tensor max_dy_tensor])
    view(-30, 30)
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$w(x,y)$')
    title(['Frame:', num2str(i), ' out of ',num2str(size(y_am,1))])
    pbaspect([a0 b0 0.5*a0])
    sgtitle('Displacement Field and Its Gradient')
    pause(0.00001)
end

%% ON DEMAND SAVING VIDEO
[ReadFileName, PrintPath] = uigetfile('*.*');
v = VideoWriter([PrintPath, 'FieldVisualization.avi']);
v.FrameRate = 10;
open(v);

figure(6),clf
for i = 1:2:size(y_am,1)
    subplot(121)
    s = surf(y_modal_tensor(:,:,i));
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    xticks(xticvec);
    yticks(yticvec);
    xticklabels(cellstr(string(xvec)));
    yticklabels(cellstr(string(yvec)));
    zlim([min_y_tensor max_y_tensor])
    xlim([1 length(x)])
    ylim([1 length(y)])
    mycolorbar();
    caxis([min_y_tensor max_y_tensor])
    view(-30, 30)
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$w(x,y)$')
    title(['Frame:', num2str(i), ' out of ',num2str(size(y_am,1))])
    pbaspect([a0 b0 0.5*a0])
    subplot(122)
    s = surf(dely_modal_tensor_tensor(:,:,i));
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    xticks(xticvec);
    yticks(yticvec);
    xticklabels(cellstr(string(xvec)));
    yticklabels(cellstr(string(yvec)));
    zlim([min_dy_tensor max_dy_tensor])
    xlim([1 length(x)])
    ylim([1 length(y)])
    mycolorbar();
    caxis([min_dy_tensor max_dy_tensor])
    view(-30, 30)
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$w(x,y)$')
    title(['Frame:', num2str(i), ' out of ',num2str(size(y_am,1))])
    pbaspect([a0 b0 0.5*a0])
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(0.00001)
end
close(v)
%% ON DEMAND SAVING GIF (FOR WEBSITE ILLUSTRATION)
PrintPath = uigetdir(cd);
filename = [PrintPath, '\FieldVisualization.gif'];
figure(6),clf
for i = 1:1:size(y_am,1)
    s = surf(y_modal_tensor(:,:,i));
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    xticks(xticvec);
    yticks(yticvec);
    xticklabels(cellstr(string(xvec)));
    yticklabels(cellstr(string(yvec)));
    zlim([min_y_tensor max_y_tensor])
    xlim([1 length(x)])
    ylim([1 length(y)])
    mycolorbar();
    caxis([min_y_tensor max_y_tensor])
    view(-30, 30)
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$w(x,y)$')
    title(['Frame:', num2str(i), ' out of ',num2str(size(y_am,1))])
    pbaspect([a0 b0 0.5*a0])
    drawnow
    frame = getframe(gcf);
    im{i} = frame2im(frame);
end

for i = 1:1:length(im)
    [imind, cm] = rgb2ind(im{i}, 256);
    if i == 1
        edit(filename)
        imwrite(imind, cm, filename,'gif','LoopCount',Inf,'DelayTime',0.0001);
    else
        imwrite(imind, cm, filename,'gif','WriteMode','append','DelayTime',0.0001);
    end
end

filename = [PrintPath, '\dFieldVisualization.gif'];
figure(7),clf
for i = 1:1:size(y_am,1)
    s = surf(dely_modal_tensor_tensor(:,:,i));
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    xticks(xticvec);
    yticks(yticvec);
    xticklabels(cellstr(string(xvec)));
    yticklabels(cellstr(string(yvec)));
    zlim([min_dy_tensor max_dy_tensor])
    xlim([1 length(x)])
    ylim([1 length(y)])
    mycolorbar();
    caxis([min_dy_tensor max_dy_tensor])
    view(-30, 30)
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$w(x,y)$')
    title(['Frame:', num2str(i), ' out of ',num2str(size(y_am,1))])
    pbaspect([a0 b0 0.5*a0])
    drawnow
    frame = getframe(gcf);
    im{i} = frame2im(frame);
end

for i = 1:1:length(im)
    [imind, cm] = rgb2ind(im{i}, 256);
    if i == 1
        edit(filename)
        imwrite(imind, cm, filename,'gif','LoopCount',Inf,'DelayTime',0.0001);
    else
        imwrite(imind, cm, filename,'gif','WriteMode','append','DelayTime',0.0001);
    end
end

%% NUMERICAL EXPERIMENT REPORT
close all
addpath('L:\My Drive\Graduate study\Research\Projects\OS')
[file, path, selectedfile] = openfiles('*.txt*');
fid = fopen(selectedfile{1},'wt');
SavePath = path;        % Pending slash
FolderName = split(SavePath,'\'); % Get the folder name and use it as the save name
SaveName = FolderName{end-1};

fprintf(fid,"====================== RESULTS FOR AMM DISCRETIZATION =================================\n");
fprintf(fid,"For Research Use: output-only modal analysis tool box v0.0\n");
fprintf(fid,"Copyright by Hewenxuan Li, hewenxuan_li@uri.edu\n");
fprintf(fid,strcat('Report Generated on: ', string(datetime),'\r\n'));
fprintf(fid,"---------------------- MEMBRANE SPECIFICATION ---------------------------------------------\n");
fprintf(fid,['Size of the membrane (a0, b0): (', num2str(a0),',', num2str(b0),')', '\n']);
fprintf(fid,['Mass Density Distribution: ', rho_xy_str, '\n']);
fprintf(fid,['Tension Distribution: ', T_xy_str, '\n']);
fprintf(fid,['Force Specification: ', F_xy_str, '\n']);
fprintf(fid,['Mass Density Constant (rho0): ', num2str(rho0), '\n']);
fprintf(fid,['Tension Constant (T0): ', num2str(T0), '\r\n']);

fprintf(fid,"-------------------------- STEPPING SETUP ---------------------------------------------\n");
fprintf(fid,['Sampling Rate (fs): ', num2str(fs), ' Hz\n']);
fprintf(fid,['Simulation Time (tend): ', num2str(tend), ' Secs\n']);
fprintf(fid,['Forcing Type: ', Forcing, '\n']);
fprintf(fid,['Stepping Method: ', Method, '\n']);
fprintf(fid,['Simulation Space: ', SimSpace, '\n']);
fprintf(fid,['Type of Excitation: ', ExciteForm, '\n']);
fprintf(fid,['Weight of Forcing: ', ForceWeight, '\r\n']);

fprintf(fid,"---------------------- DISCRETIZATION INFORMATION -------------------------------------\r\n");
fprintf(fid,['Indices of the x-direction eigenfunctions considered (m): ', num2str(1:mm), '\n']);
fprintf(fid,['Indices of the y-direction eigenfunctions considered (n): ', num2str(1:nn), '\n']);
fprintf(fid,['Number of Modes Considered (nm): ', num2str(nm), ' \n']);
fprintf(fid,['Spatial Resolution (dx, dy): (', num2str(dx),', ', num2str(dy),')']);

fprintf(fid,'The stiffness matrix (first 5-by-5):\n');
for k = 1:5
    fprintf(fid,' |%.2f         %10.2f          %10.2f          %10.2f          %10.2f|\n', K_am(k,1:5));
end
PathStr = split(path,'\');
PathStr = strjoin(PathStr,'//');
fprintf(fid,'\n');
fprintf(fid,['See full systems matrix in:', PathStr, 'Data.m\n']);
fprintf(fid,'using Data.M, Data.C, and Data.K\n');
fprintf(fid,'\n');
fprintf(fid,'The stiffness matrix (first 5-by-5):\n');
for k = 1:5
    fprintf(fid,' |%.2f         %10.2f          %10.2f          %10.2f          %10.2f|\n', M_am(k,1:5));
end
fprintf(fid,'\n');
fprintf(fid,['Type of Damping:', DampingType, '\n']);
fprintf(fid,'Modal Damping (mode 1-5): |%.3f, %.3f, %.3f, %.3f, %.3f|\n', zeta_am(1:5));
fprintf(fid,'The damping matrix (first 5-by-5):\n');
for k = 1:5
    fprintf(fid,' |%.2f         %10.2f          %10.2f          %10.2f          %10.2f|\n', C_am(k,1:5));
end
fprintf(fid,'\n');
fprintf(fid,'The forcing weight (first 5 modes):\n');
fprintf(fid,' |%.2f         %10.2f          %10.2f          %10.2f          %10.2f|\n', weight(1:5));
fprintf(fid,'\n');
fprintf(fid,'The natural frequencies (rad/sec):\n');
fprintf(fid,' |%.3f         %10.3f          %10.3f          %10.3f          %10.3f|\n', omega_am(1:5));
fprintf(fid,'\n');
fprintf(fid,'The damped natural frequencies (rad/sec):\n');
fprintf(fid,' |%.3f         %10.3f          %10.3f          %10.3f          %10.3f|\n', omegad_am(1:5));
fprintf(fid,'\n');

fprintf(fid,"---------------------- RESPONSE INFORMATION -------------------------------------\r\n");
fprintf(fid,'\n');
fprintf(fid,['The Size of the Configuration Data:',num2str(size(y_modal,1)),'-by-',num2str(num2str(size(y_modal,2))),'\n']);
fprintf(fid,['The Manitude of Modal Participation (first 5 modes):',num2str(IW(2,1:5)),'\n']);
fprintf(fid,['Modes Included: ', num2str(ModeIndices), '\n']);
fprintf(fid,"=======================================================================================");

fclose(fid);
disp('Report generated!')
%% SAVE THE DATA INTO A DATA STRUCTURE
Data.System = 'Non-uniform Fixed-Fixed Membrane'; % System Type
Data.Date = string(datetime);                  % Create Time Stamp
Data.DiscretizationMethod = 'AMM';             % Discretization Method Used
% System Specification
Data.Dimensions.a0 = a0;                       % Membrane Dimension in x
Data.Dimensions.b0 = b0;                       % Membrane Dimension in y
Data.SpatialResolution.dx = dx;                % Spatial Resolution in x
Data.SpatialResolution.dy = dy;                % Spatial Resolution in y
Data.Density = rho_xy_str;                     % Nonuniform Density (definition)
Data.Tension = T_xy_str;                       % Nonuniform Tension (definition)
Data.Force = F_xy_str;                         % Forcing (definition)
Data.UniformDensity = rho0;                    % Uniform Material Density
Data.UniformTension = T0;                      % Uniform Tension
% Stepping Specification
Data.fs = fs;                                  % Sampling Frequency
Data.tend = tend;                              % Ending Time for Simulation
Data.Forcing = forcing;                        % Forcing Function
Data.SimulationMethod = Method;                % Simulation Method
Data.SimulationSpace = SimSpace;               % Simulation Space
Data.ExciteForm = ExciteForm;                  % Excitation Form
Data.ForceWeight = ForceWeight;                % Force Weight
Data.IW = IW;                                  % Initial Modal Participation
Data.ModeIndices = ModeIndices;                % Mode Indices
% Discretization Setup
Data.xEigFcn = 1:mm;                           % Assign Eigenfunctions in x direction
Data.yEigFcn = 1:nn;                           % Assign Eigenfunctions in y direction
Data.nm = nm;                                  % Number of modes in the simulation
% Discretization Results
Data.M = M_am;                                 % Mass matrix of the discretized system
Data.K = K_am;                                 % Stiffness matrix of the discretized system
Data.C = C_am;                                 % Damping matrix of the discretized system
Data.weight = weight;                          % Weight to the applied force for the discretized system
Data.omegan = omega_am;                        % Natural frequencies
Data.omegad = omegad_am;                       % Damped natural frequencies 
% Mode Shapes  
Data.Phi_x = Phi_x;                            % Displacement mode shapes 
Data.Phi_v = Phi_v;                            % Velocity mode shapes 
% Response Data
Data.Eta = eta_am;                             % Simulated modal motion
Data.Y = y_modal;                              % Simulated configuration motion
Data.delY = dely_modal;                        % Simulated gradient motion
Data.TY = y_modal_tensor;                      % Simulated configuration motion tensor
Data.TdelY = dely_modal_tensor;                % Simulated modal motion tensor
% Response Size
Data.SnapSize = [size(y_modal_tensor,1), size(y_modal_tensor,1)]; % Size of each of the snapshot of the field
Data.ResponseMatSize = size(y_modal);          % Size of the response matrix
% Save Data to Data.mat
save([path,'Data.mat'],'Data')
disp(strcat('Data Saved! ',string(datetime)))
%%
x = 0:0.01:1;
n = 3;
y = x.^n;
plot(x,y,'LineWidth',3)
title(['y = x^n,  n = ' num2str(n) ])

n = 1:0.5:5;
nImages = length(n);

fig = figure;
for idx = 1:nImages
    y = x.^n(idx);
    plot(x,y,'LineWidth',3)
    title(['y = x^n,  n = ' num2str( n(idx)) ])
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end
close;

figure;
for idx = 1:nImages
    subplot(3,3,idx)
    imshow(im{idx});
end

% filename = 'testAnimated.gif'; % Specify the output file name
for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end
%% FOR NOW, THIS CODE ONLY WORKS FOR LINEAR SYSTEMS!
% Prepare for the substitution of the derivative equation of motion
dF = ppval(fnder(forcing),tspan);
dF = repmat(dF,[10 1]);
W = repmat(weight, [1 size(dF,2)]);
dF = W.*dF;
% Back plug the state vector y to obtain the velocity state vector
dq_am = A_am*q_am' + B*Um';
dq_am = dq_am'; % derivative data
% Back substitute the velocity into the derivative equation of motion
ddq_am = A_am*dq_am' + B*dF;
ddq_am = ddq_am'; % 

% Rearrange the response data such that only one component remains
ddq_am = dq_am(:,nm+1:2*nm);  % Acceleration
dddq_am = dq_am(:,nm+1:2*nm); % Jerk
q_am = q_am(:,1:nm);          % Disp
dq_am = q_am(:,nm+1:2*nm);    % Velocity

ddy_am = ddq_am*Phi_am';
dddy_am = dddq_am*Phi_am';

%% Animation of the modal motion
close all
eta = real(eta_am(:,1:1:end));
disp('MODAL PARTICIPATION OVER THE DESIGNATED TIME FRAME:')
f1 = figure;
f1.Renderer = 'opengl';
str = 1;
pr = 20;                             % Playback range (time stamp)
dpr = 1;                        % Playback resolution (time resolution)
autorange = max(max(eta(str:dpr:str+pr,1:n)));                  % Auto-ranging the plot range
clf
for i = str:dpr:str+pr
    stem([eta(i,1:n)],'Color','r')                  % Plot real-time modal participation
    pause(0.005)
    
    hold on
    stem([eta(i,1:n)],'Color','#D0D0D0')            % Hold on the trace of the modal participation
    ylim([-autorange*1.5 autorange*1.5]) 
    if i == 1
        xlabel('Mode Index')
        ylabel('Amplitude')
        title('Snapshot of the Modal Participation')
        pbaspect([1.618 1 1])
%         xticks([1 2 3 4 5 6 7 8 9 10])
%         xticklabels({'#1','#2','#3','#4','#5','#6','#7','#8','#9', '#10'})
    end
end
disp('Modal Participation Visualization Done!')
stem([eta(i,1:n)],'r')                              % Retain the last played position

%% PREVIEW: Animation of the configuration motion
clf
% axis tight manual 
set(gca,'nextplot','replacechildren'); 
set(gcf,'renderer','painters')
str = 1;
pr = 200;                                                  % Playback range (time stamp)
dpr = 1;                                                  % Playback resolution (time resolution)
plength = length(str:dpr:str+pr);                         % Playback time vector length
PlotData = y_am;
autorange = max(max(PlotData(str:dpr:str+pr,1:size(y_am,2))));           % Auto-ranging the plot range
c = linspace(0,0.6,plength);     
% Color gradient
init_color = [0.9 0.9 0.9];                               % Initial color
counter = 0;
clf
for i = str:dpr:str+pr
    counter = counter + 1;
    plot([0 PlotData(i,:)], 'k', 'LineWidth', 2)                           % Plot real-time configuration space motion (zero padded)
    xlim([0 length(x) + length(x)*0.05])
    ylim([-autorange*1.2 autorange*1.2]) 
    pbaspect([2 1 1])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(['Frame: ',num2str(i),' out of ',num2str(str+pr),'.'])
%     axis off
    set(gcf, 'papersize', [6 6*0.5])
    set(gcf, 'paperposition', [0 0 6 6*0.5])
    pause(dt*dpr)
end
%% Animation of the configuration motion
% clf
% SavePath = uigetdir('','Set the path to which the trimmed figure is saved.');
% SavePath = [SavePath,'\'];
% SaveName = 'Beam_Harmonic_Mode_1.avi';
% % axis tight manual 
% set(gca,'nextplot','replacechildren'); 
% set(gcf,'renderer','painters');
% % Video Writer
% v = VideoWriter([SavePath, SaveName],'Uncompressed AVI');
% v.FrameRate = 10;
% open(v)
% str = 1;
% pr = 200;                                                % Playback range (time stamp)
% dpr = 1;                                                  % Playback resolution (time resolution)
% plength = length(str:dpr:str+pr);                         % Playback time vector length
% PlotData = y_am;
% autorange = max(max(PlotData(str:str+pr,1:size(y_am,2))));           % Auto-ranging the plot range
% c = linspace(0,0.6,plength);     
% % Color gradient
% init_color = [0.9 0.9 0.9];                               % Initial color
% counter = 0;
% clf
% for i = str:dpr:str+pr
%     counter = counter + 1;
%     plot([0 PlotData(i,:)], 'k', 'LineWidth', 2)                           % Plot real-time configuration space motion (zero padded)
%     xlim([0 length(x) + 1])
%     ylim([-autorange*1.1 autorange*1.1]) 
%     pbaspect([2 1 1])
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     set(gcf, 'papersize', [6 6*0.5])
%     set(gcf, 'paperposition', [0 0 6 6*0.5])
%     title([SaveName,'; frame: ',num2str(i),' out of ',num2str(str+pr),'.'],'Interpreter','latex')
%     frame = getframe(gcf);
%     writeVideo(v,frame);
%     pause(dt*dpr)
% end
% close(v);
%% Save Animated Pictures 
clf
SavePath = uigetdir('','Set the path to which the trimmed figure is saved.');
SavePath = [SavePath,'\'];
SaveName = 'Beam_Harmonic_Thick_Mode_2&3';
ExtName = '.png';
% axis tight manual 
set(gca,'nextplot','replacechildren'); 
set(gcf,'renderer','painters');
% Video Writer
str = 1;
pr = 2000;                                                % Playback range (time stamp)
dpr = 1;                                                  % Playback resolution (time resolution)
plength = length(str:dpr:str+pr);                         % Playback time vector length
PlotData = y_am;
autorange = max(max(PlotData(str:str+pr,1:size(y_am,2))));           % Auto-ranging the plot range
c = linspace(0,0.6,plength);     
% Color gradient
init_color = [0.9 0.9 0.9];                               % Initial color
counter = 0;
tic
for i = str:dpr:str+pr
    counter = counter + 1;
    plot([0 PlotData(i,:)], 'k', 'linewidth', 2)                           % Plot real-time configuration space motion (zero padded)
    xlim([0 length(x) + length(x)*0.05])
    ylim([-autorange*1.2 autorange*1.2])
    pbaspect([2 1 1])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gcf, 'papersize', [6 6*0.5])
    set(gcf, 'paperposition', [0 0 6 6*0.5])
    title([SaveName,'; frame: ',num2str(i),' out of ',num2str(str+pr),'.'],'Interpreter','latex')
    pause(dt*dpr)
    print([SavePath,SaveName,'_',num2str(i),ExtName],'-dpng','-r150')
end
toc
%% Animation of the configuration motion
set(gcf,'renderer','opengl')
str = 1;
pr = 100;                             % Playback range (time stamp)
dpr = 1;                        % Playback resolution (time resolution)
plength = length(str:dpr:str+pr);                         % Playback time vector length
PlotData = y_am;
autorange = max(max(PlotData(str:str+pr,1:size(y_am,2))));                  % Auto-ranging the plot range
c = linspace(0,0.6,plength);     
% Color gradient
init_color = [0.9 0.9 0.9];                            % Initial color
counter = 0;
clf
for i = str:dpr:str+pr
    counter = counter + 1;
    plot([0 PlotData(i,:)], 'r')                           % Plot real-time configuration space motion (zero padded)
    pause(0.0025)
    hold on
    plot([0 PlotData(i,:)],'linestyle','-','color',(init_color)*(1-c(counter)))  % Persistent plot of the trace of the motion
    if counter == 1
        xlabel('Nodal Position')
        ylabel('Amplitude')
        title('Snapshot of the Beam Vibration')
        pbaspect([1.618 1 1])
        axis tight
%         xticks([1 2 3 4 5 6 7 8 9 10 11])
%         xticklabels({'Base','#1','#2','#3','#4','#5','#6','#7','#8','#9', '#10'})
    end
    ylim([-autorange*1.5 autorange*1.5]) 

end
plot([0 PlotData(i,:)], 'r','linewidth',1)

%% Plot of the modal participation through their Standard Deviations
figure(5),clf
str = 1;                             % Playback range (time stamp)
stp = 23140;                             % Playback range (time stamp)
autorange = max((std(q_am(:,1:n))));                  % Auto-ranging the plot range
stem(std(q_am(str:stp,1:n)),'Color','r')                  % Plot real-time modal participation
ylim([0 autorange*1.5])
xlabel('Index of Mode')
ylabel('std(Modal Response)')
title('Snapshot of the Modal Participation')
pbaspect([1.618 1 1])

% FRF Estimation
disp('THE FOLLOWING FIGURES ARE BELONGS TO THE FRF ESTIMATION SECTION!')
figure(6),clf
str = 1;                             % Playback range (time stamp)
stp = 15250;                             % Playback range (time stamp)
FreqRangeMin = 0;
FreqRangeMax = 155;                        % Frequency display control (xlim)
YScaling = "log";
Channels = 4;
windsize = 4096;
noverlap = windsize/2;
signal = q_am;
% Calculation: FRF and Pyy
[FRF, f] = modalfrf(ppval(forcing,tspan(str:stp)), signal(str:stp,Channels)-mean(signal(str:stp,Channels)), fs, boxcar(windsize), noverlap, 'Sensor', 'dis');
[pyy, ~] = pwelch(signal(str:stp,Channels)-mean(signal(str:stp,Channels)), boxcar(windsize), noverlap, [], fs);

% --------------------- VISUALIZATION BEGINS FROM HERE ------------------
% TIME SERIES AND PSD FIGURES BEGIN HERE
figure(10), clf
subplot(211)
plot(signal(str:stp,Channels))
axis tight
xlabel('Sample')
ylabel('Magnitude')
title('Truncated Signal') 
pbaspect([3 1 1])
subplot(212)

if isequal(lower(YScaling), 'log')
    plot(f, 10*log10(pyy))
    hold on
    for i = 1:length(omegad_am)
        plot([omegad_am(i) omegad_am(i)]/2/pi, [min(min(10*log10(pyy))) max(max(10*log10(pyy)))], 'k:') % MODE LINES
    end
    if isequal(lower(Forcing),'harmonic')
        plot([finfo.omega finfo.omega]/2/pi, [min(min(10*log10(pyy))) max(max(10*log10(pyy)))], 'r--')  % FORCING LINE
    end
else
    plot(f, pyy)
    hold on
    for i = 1:length(omegad_am)
        plot([omegad_am(i) omegad_am(i)]/2/pi, [min(min(pyy)) max(max(pyy))], 'k:') % MODE LINES
    end
    if isequal(lower(Forcing),'harmonic')
        plot([finfo.omega finfo.omega]/2/pi, [min(min(pyy)) max(max(pyy))], 'r--')  % FORCING LINE
    end
end

hold off

xlabel('Frequency (Hz)')
ylabel('PSD')
title('Power Spectral Density Estimate')
pbaspect([3 1 1]), xlim([FreqRangeMin FreqRangeMax]), 

if isequal(lower(YScaling), 'log')
    ylim([min(min(10*log10(pyy))) max(max(10*log10(pyy)))]); % Plot setting
else
    ylim([min(min(pyy)) max(max(pyy))]); % Plot setting
end

% FREQUENCY RESPONSE FUNCTION BEGINS THERE IF IT IS NOT FREEDECAY
if ~isequal(lower(Forcing),'freedecay')
    
% FRF FIGURES BEGIN HERE
figure(11), clf
subplot(211)
for i = 1:length(Channels)
    if isequal(lower(YScaling), 'log')
        if i<=4
            pfrf = plot(f,log10(abs(FRF(:,i))), 'LineStyle',LineStyles{i}, 'Color', Mycolors.defaultcolor{1});
        elseif i > 4 && i <= 8
            pfrf = plot(f,log10(abs(FRF(:,i))), 'LineStyle',LineStyles{i-4}, 'Color', Mycolors.defaultcolor{1});
            hold on
            plot(f(1:10:end), log10(abs(FRF(1:10:end,i))),Markers{i-4},'MarkerSize',3, 'MarkerEdgeColor', Mycolors.defaultcolor{1})
        elseif i > 8
            pfrf = plot(f,log10(abs(FRF(:,i))), 'LineStyle',LineStyles{i-8}, 'Color', Mycolors.defaultcolor{1});
            hold on
            plot(f(1:10:end), log10(abs(FRF(1:10:end,i))),Markers{i-4},'MarkerSize',3, 'MarkerEdgeColor', Mycolors.defaultcolor{1})
        end
    else
        if i<=4
            pfrf = plot(f,abs(FRF(:,i)), 'LineStyle',LineStyles{i}, 'Color', Mycolors.defaultcolor{1});
        elseif i > 4 && i <= 8
            pfrf = plot(f,abs(FRF(:,i)), 'LineStyle',LineStyles{i-4}, 'Color', Mycolors.defaultcolor{1});
            hold on
            plot(f(1:10:end), abs(FRF(1:10:end,i)),Markers{i-4},'MarkerSize',3, 'MarkerEdgeColor', Mycolors.defaultcolor{1})
        elseif i > 8
            pfrf = plot(f,abs(FRF(:,i)), 'LineStyle',LineStyles{i-8}, 'Color', Mycolors.defaultcolor{1});
            hold on
            plot(f(1:10:end), abs(FRF(1:10:end,i)),Markers{i-4},'MarkerSize',3, 'MarkerEdgeColor', Mycolors.defaultcolor{1})
        end
%         pfrf = plot(f,abs(FRF(:,i)), 'LineStyle',"-",'Marker',Markers{i},'MarkerSize',3);
    end
hold on
end

if isequal(lower(YScaling), 'log')
    for i = 1:length(omegad_am)
        ml = plot([omegad_am(i)/2/pi omegad_am(i)/2/pi], [min(min(log10(abs(FRF)))) max(max(log10(abs(FRF))))], 'k:'); % MODE LINES
    end
    if isequal(lower(Forcing),'harmonic')
        fl = plot([finfo.omega finfo.omega]/2/pi, [min(min(log10(abs(FRF)))) max(max(log10(abs(FRF))))], 'r--');           % FORCING LINE
    end
else
    for i = 1:length(omegad_am)
        ml = plot([omegad_am(i)/2/pi omegad_am(i)/2/pi], [min(min(abs(FRF))) max(max(abs(FRF)))], 'k:'); % MODE LINES
    end
    if isequal(lower(Forcing),'harmonic')
        fl = plot([finfo.omega finfo.omega]/2/pi, [min(min(abs(FRF))) max(max(abs(FRF)))], 'r--');           % FORCING LINE
    end
end
xlabel('Frequency (Hz)'), ylabel('FRF'), title('Frequency Response Function Estimate') 
pbaspect([3 1 1])
xlim([FreqRangeMin FreqRangeMax])
if isequal(lower(YScaling), 'log')
    ylim([min(min(log10(abs(FRF)))) max(max(log10(abs(FRF))))]); % Plot setting
else
    ylim([min(min(abs(FRF))) max(max(abs(FRF)))]); % Plot setting
end
if isequal(lower(Forcing),'harmonic')
    legend([pfrf, ml, fl],{'FRF','$f_n$','$f_f$'})
else
    legend([pfrf, ml],{'FRF','$f_n$'})
end
subplot(212)
% Plot individual FRF
for i = 1:size(FRF,2)
    if i<=4
        php = plot(f,wrapToPi(angle(FRF(:,i))), 'LineStyle', LineStyles{i}, 'Color', Mycolors.defaultcolor{2});
        hold on
    elseif i>4 && i<=8
        plot(f,wrapToPi(angle(FRF(:,i))), 'LineStyle', LineStyles{i-4}, 'Color', Mycolors.defaultcolor{2});
        hold on
        plot(f(1:10:end), wrapToPi(angle(FRF(1:10:end, i))),Markers{1},'MarkerSize',3, 'MarkerEdgeColor', Mycolors.defaultcolor{2})
    elseif i > 8
        plot(f,wrapToPi(angle(FRF(:,i))), 'LineStyle', LineStyles{i-8}, 'Color', Mycolors.defaultcolor{2});
        hold on
        plot(f(1:10:end), wrapToPi(angle(FRF(1:10:end, i))),Markers{1},'MarkerSize',3, 'MarkerEdgeColor', Mycolors.defaultcolor{2})
    end
end

for i = 1:length(omegad_am)
    ml = plot([omegad_am(i)/2/pi omegad_am(i)/2/pi], [min(min(wrapToPi(angle(FRF)))) max(max(wrapToPi(angle(FRF))))], 'k:'); % MODE LINES
end
if isequal(lower(Forcing),'harmonic') 
    fl = plot([finfo.omega finfo.omega]/2/pi, [min(min(wrapToPi(angle(FRF)))) max(max(wrapToPi(angle(FRF))))], 'r--');           % FORCING LINE
end
title('Phase of the Estimated FRF')
xlabel('Frequency (Hz)')
ylabel('Phase Angle (rad)')
set(gca,'color','none')
xlim([FreqRangeMin FreqRangeMax])
ylim([min(min(wrapToPi(angle(FRF)))) max(max(wrapToPi(angle(FRF))))])

end

Plot the Assumed Modes
ModeIndex = [3, 5, 7, 9];
disp('THE FOLLOWING FIGURE BELONGS TO THE MODE SHAPE SECTION:')
figure(20), clf
plot(phi_am(:, ModeIndex)./phi_am(end, ModeIndex))
hold on
plot([0 11],[0 0], 'k:')
xlim([1 11])
ylim([-1.1 1.1])
pbaspect([2 1 1])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Base','#1','#2','#3','#4','#5','#6','#7','#8','#9', '#10'})
grid on
xlabel('Position index of the Discretized Beam')
ylabel('Tip-Deflection-Normalized Magnitude')

% Save the data
% Select 'True' from the drop-down bar 
% A pop up window is shown, during which the users have to designate the location to save the file
% 
% This saving function follows a naming rule:
% The file name starts with 'CantBeam_simulation_'
% If the forcing is applied from the base of the system: 
% 'CantBeam_simulation' + 'Forcing Type' + 'Forcing Information' + 'Sampling Rate' + 'Simulation Duration' + 'Force Weight (weighted by all modes or a single point on the beam)' + 'Simulation Method'
% As the forcing type veries, the 'Forcing Information' changes: if it is random, then this information include only the std of the random loading; if the loading is harmonic, the information include both amplitude and frequency; if it is freedecay, the information is empty.
% When the freedecay is selected, the 'Force Weght' is also empty; instead, one expects which mode is present and how strong the modes are participated.
% If the forcing is applied from the beam (at some arbitrary location): 
% 'CantBeam_simulation' + 'Forcing Type' + 'Forcing Information' + 'Sampling Rate' + 'Simulation Duration' + 'Force Weight (weighted by all modes or a single point on the beam)' + 'Excitation Location' + 'Simulation Method'; the 'excitation location' is added to the saving name.
% 
% Sample Names:
%     CantBeam_simulation_FreeDecay_1000_200_1  2  3  4_1_LSIM.mat
% Indicates the simulation is of 'Freedecay' type + the sampling frequency is 1000Hz + simulation duration is 200 seconds + there are four modes got involved, namely, [1 2 3 4], with initial weight 1 (each mode has identical contribution) + simulation method 'LSIM' + proportionally damped.
%     CantBeam_simulation_Impulse_1000_200_10_2_LSIM_0.05_NPD.mat
% Indicates the simulation of 'Impulse' response + the sampling rate is 1000Hz + simulation duration is 200 seconds + the excitation point is located at the 10th node + 2 impulses in total + using the LSIM Function to simulate the system + 0.05 externally attached damping + nonproportional damping.

SaveData = 'True';
% Save name convention: 
% CantBeam_simulation_ForcingType_ForcingAmp_SamplingFreq_TotalTime_ForcingLocation(Beam/Base)_ForcingWeight_SimulationMethod
if isequal('true', lower(SaveData))
    if isequal(lower(ForceWeight), 'base')
        if isequal(lower(Forcing), 'random')
            SaveName = ['CantBeam_simulation_',Forcing,...
                '_',num2str(finfo.A),'_',num2str(fs),'_',num2str(timing.Tend),...
                '_',num2str(ForceWeight),'_',num2str(Method)];
        elseif isequal(lower(Forcing), 'harmonic')
            SaveName = ['CantBeam_simulation_',Forcing,...
                '_',num2str(finfo.A),'_',num2str(finfo.omega),'_',num2str(fs),...
                '_',num2str(timing.Tend),'_',num2str(ForceWeight),'_',num2str(Method)];
        elseif isequal(lower(Forcing),'freedecay')
            SaveName = ['CantBeam_simulation_',Forcing,'_',num2str(fs),...
                '_',num2str(timing.Tend),'_',num2str(q0indx),'_',...
                num2str(InitialWeight),'_',num2str(Method)];
        elseif isequal(lower(Forcing),'impulse')
            SaveName = ['CantBeam_simulation_',Forcing,'_',num2str(fs),...
                '_',num2str(timing.Tend),'_',num2str(ExciteLoc),'_',num2str(finfo.Nimp),'_',num2str(Method)];
        elseif isequal(lower(Forcing), 'burstrandom')
                SaveName = ['CantBeam_simulation_',Forcing,...
                    '_',num2str(finfo.A),'_',num2str(fs),'_',num2str(timing.Tend),...
                    '_',num2str(ForceWeight),'_',num2str(finfo.N/2),'_',num2str(Method)];
        end
    else
        if isequal(lower(Forcing), 'random')
            SaveName = ['CantBeam_simulation_',Forcing,'_',num2str(finfo.A),...
                '_',num2str(fs),'_',num2str(timing.Tend),'_',...
                num2str(ForceWeight),'_',num2str(ExciteLoc),'_',num2str(Method)];
        elseif isequal(lower(Forcing), 'harmonic')
            SaveName = ['CantBeam_simulation_',Forcing,...
                '_',num2str(finfo.A),'_',num2str(finfo.omega),'_',num2str(fs),...
                '_',num2str(timing.Tend),'_',num2str(ForceWeight),'_',...
                num2str(ExciteLoc),'_',num2str(Method)];
        elseif isequal(lower(Forcing),'freedecay')
            SaveName = ['CantBeam_simulation_',Forcing,'_',num2str(fs),...
                '_',num2str(timing.Tend),'_',num2str(q0indx),'_',...
                num2str(InitialWeight),'_',num2str(Method)];
        elseif isequal(lower(Forcing),'impulse')
            SaveName = ['CantBeam_simulation_',Forcing,'_',num2str(fs),...
                '_',num2str(timing.Tend),'_',num2str(ExciteLoc),'_',num2str(finfo.Nimp),'_',num2str(Method)];
        elseif isequal(lower(Forcing), 'burstrandom')
            SaveName = ['CantBeam_simulation_',Forcing,...
                '_',num2str(finfo.A),'_',num2str(fs),'_',num2str(timing.Tend),...
                '_',num2str(ForceWeight),'_',num2str(ExciteLoc),'_',num2str(finfo.N/2),'_',num2str(Method)];
        end
    end
    
    if zeta_am ~= 0
        SaveName = [SaveName, '_zeta_', num2str(zeta_am(1))]
    end
    
    if M0 ~= 0
       SaveName = [SaveName,'_M0_',num2str(M0)]; 
    end
    
    if k0 ~= 0
        SaveName = [SaveName,'_k0_',num2str(k0)];
    end
    
    if npd == 1
        SaveName = [SaveName, '_C_', num2str(C0),'_iC0_',num2str(iC0),'_NPD'];
    end
    
    SavePath = uigetdir('',['Set the path to which the vibration data is saved. File name:', SaveName]);
    SavePath = [SavePath,'\'];
    % Create the saving package:
    simulation.response = {y_am, dy_am, ddy_am, dddy_am};    % Assign response to the simulation pack
    simulation.modal_response = {q_am, dq_am, ddq_am, dddq_am}; % Assign modal response to the pack
    simulation.forcing = forcing;  % Assign forcing to the simulation pack
    simulation.finfo = finfo;      % Assign finfo info to the sim. pack
    simulation.timing = timing;    % Assign timing info to the sim. pack
    simulation.system = {A_am, Km_am, Cm_am, Phi_am, zeta_am}; % Assign system info to the sim. pack
    if isequal(lower(Forcing),'freedecay')
        simulation.q0indx = q0indx;
        simulation.q0 = q0;
    end
    save([SavePath,SaveName,'.mat'],'simulation')
end
disp([savepath,SaveName,'.mat'])
