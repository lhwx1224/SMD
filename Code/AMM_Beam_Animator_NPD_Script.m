% Forced Beam Assumed Modes Method (AMM) - Linear cantilever beam simulation
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
addpath('Output-only Modal Analysis Toolbox','OS')
[Markers, LineStyles, Mycolors] = mystyle();   % Loading printing and plotting style
Forcing = 'FreeDecay';              % External Excitation Form
Method = 'LSIM';               % Simulation Method
ExciteForm = 'Beam';     % Excitation Formulation (Beam - original, Base - d^2(original)/dt^2)
ForceWeight = 'Beam';    % Weight used in weighting the force in the AMM
ExciteLoc = 10; % Excitation location
% Initial value problem setup
IW = zeros(10,1);                      % Initial weight vector
InitialWeight = 1;                     % Initial weight value
q0indx = 1:10;         % To which mode the weight is assigned
IW(q0indx) = InitialWeight;
IW = [1; 0.8; 0.6; 0.125; 0.065; 0.03; 0.015; 0.0075; 0.003; 0.001];
% Simulation setup
L = 1;                              % Length of the beam
m = 1;                              % Mass density 
M0 = 0;
k0 = 0;
EI = 1;                             % Flexural rigidity
dx = 0.001*L;                         % sampling spatial distance
x = 0:dx:L;                         % define domain of the beam
dt = 0.01;                          % sampling time
fs = 1/dt;                          % sampling frequency
tend = 5;                         % ending time of simulation
t = 0:dt:tend;                      % Time vector
nm = 5;                             % Number of modes to consider
rr = 2;                             % set the resapling rate (optional). Test other than random should use rr = 1;

% Eigenvalues to the DEVP (From Mathematica Calculation)
beta = [1.8751, 4.69409, 7.85476, 10.9955, 14.1372, ...
    17.2788, 20.4204, 23.5619, 26.7035, 29.8451, ...
    32.9867, 36.1283, 39.2699, 42.4115, 45.5531, ...
    48.6947, 51.8363, 54.9779, 58.1195, 64.4026]';
% Amplitude to normalize the modes (From Mathematica Calculation)
A = [1.36222, 0.981868, 1.00078, 0.999966, 1, ...
    1, 1, 1., 0.999999, 1, ...
    0.999853, 1.00246, 0.993615, 0.960281, 0.920944, ...
    11952.2, 0.864524, 0.839219, 0.815037, 0.775349]';
% Modal Damping Ratios
zeta_am = 0.00*ones(nm,1);
zeta_am = 0.001*(1:1:10);

a = (sin(beta(1:nm))+sinh(beta(1:nm)))./(cos(beta(1:nm)) + cosh(beta(1:nm)));
omega_am = beta(1:nm).^2;
gamma_am = zeta_am.*omega_am;
omegad_am = omega_am.*sqrt(1 - zeta_am.^2);
phi_am = (sin(beta(1:nm)*x) - sinh(beta(1:nm)*x) - diag(a)*(cos(beta(1:nm)*x) - cosh(beta(1:nm)*x)))';
% normalize modes to unit norm
phi_am = phi_am/diag(A(1:nm));
Phi_am = phi_am(2:end,:);
% Build Matrices (Dynamic matrices in configuration space)
n = nm;                                % number of nodes
Mm_am = eye(n);                        % Mass matrix (eye if uniform beam with cant beam assumption)
Mm_am = Mm_am + M0*Phi_am(n, :)'*Phi_am(n,:);
Km_am = diag(omega_am.^2);             % Stiffness matrix
Km_am = Km_am + k0*Phi_am(n, :)'*Phi_am(n,:);
Cm_am = diag(2*gamma_am);              % Damping matrix
% Add a damper in the middle of the beam
C0 = 0; % damping coefficient of the damper in the middle
iC0 = 5; % location of the lumped damper
% Append non-proportional damping effect to the damping matrix
Cm_am = Cm_am + C0*Phi_am(iC0, :)'*Phi_am(iC0,:);
% Natural modes from the mass and stiffness matrices
[U, lambda] = eig(Mm_am\Km_am);
% State-space Dynamic Matrix
A_am = [zeros(n), Mm_am; -Km_am, -Cm_am]; 

% Check if the system is proportionally damped:
if sum(Km_am*eye(n)*Cm_am == Cm_am*eye(n)*Km_am,"all") == n^2
    disp('System is proportionally damped!')
    npd = 0;
else
    disp('System is non-proportionally damped!')
    npd = 1;
end
% ===================== Get forcing function =============================
switch lower(Forcing)
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
        DataInspection(ppval(forcing,tspan),tspan,fs)
        pause(0.1)
end

disp(['ADJUSTED SAMPLING FREQUENCY:', num2str(fs),'Hz'])

% Weight in front of the forcing vector assuming base motion (from Mathematica)
if isequal(lower(ForceWeight), "base")
    weight = [0.782992, 0.433936, 0.254425, 0.181898, 0.141471, 0.115749, ...
        0.0979415, 0.0848826, 0.0748957, 0.0670128]';
else
    weight = Phi_am(ExciteLoc,:)';
end

B = [zeros(n,n);eye(nm)];   % Forcing weight matrix

if isequal(lower(Forcing),'freedecay')
    q0 = zeros(2*n,1);            % Initial conditions
    q0(1:n,1) = IW(1:n)+q0(1:n,1);
else
    q0 = zeros(2*n,1);    % Initial conditions (default zeros(2*n,1))   
end
% q0 = [sum(Phi_am(:,1),2); zeros(n,1)];            % Initial conditions
% q0 = zeros(2*n,1);            % Initial conditions


switch lower(Method)
    case 'lsim'
 % ------------------- Linear System Simulation "lsim" -------------------
        Fm = repmat(ppval(forcing,tspan),nm,1)'; % Augmented forcing matrix
        Um = zeros(size(Fm));                    % allocate weighted forcing matrix
        for i = 1:length(weight) % The for loop can be replaced by a repmat of weight
            Um(:,i) = weight(i)*Fm(:,i); % weighted augmented forcing matrix
        end
        sys = ss(A_am, B, eye(size(A_am)), zeros(size(B))); % Define state space system in transfer function form
        tic
        [q_am, t] = lsim(sys, Um, tspan, q0);                      % Solve the system
        toc
        disp('INSPECTION OF THE MODAL COORDINATES:')
        DataInspection(q_am(:,1:nm),t,fs)
        % y_am = q_am(:,1:n)*U'*Phi_am';
        y_am = q_am(:,1:n)*Phi_am';
        dy_am = q_am(:,n+1:2*n)*Phi_am';
        disp('INSPECTION OF THE PHYSICAL COORDINATES:')
        DataInspection(y_am,t,fs)
% -------------------- Simulation via Numerical Integration --------------
    case 'ode'
        % Numerical integration
        odefcn = @(t,q) GenLinBeam_conf(t,q,A_am,B,forcing,weight); % DEFINE ODE FUNCTION
        % options = odeset('OutputFcn',@(t,y,flag) odeoutput_bw(t,y,flag,m,a,k,B,alpha,beta1,beta2,gamma,c,n,f));
        disp('Numerical integration started:')
        tic
        [t,q_am] = ode45(odefcn,tspan,q0);
        toc
        disp('Numerical integration finished!')
        % y_am = q_am(:,1:n)*Phi_am';
        dy_am = q_am(:,n+1:2*n)*Phi_am';
        y_am = q_am(:,1:n)*Phi_am';
        DataInspection(q_am(:,1:nm),t,fs)
        DataInspection(y_am,t,fs)
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
disp('MODAL PARTICIPATION OVER THE DESIGNATED TIME FRAME:')
str = 1;
pr = 2000;                             % Playback range (time stamp)
dpr = 17;                        % Playback resolution (time resolution)
autorange = max(max(q_am(str:dpr:str+pr,1:n)));                  % Auto-ranging the plot range
clf
for i = str:dpr:str+pr
    stem([0 q_am(i,1:n)],'Color','r')                  % Plot real-time modal participation
    pause(0.005)
    
    hold on
    stem([0 q_am(i,1:n)],'Color','#D0D0D0')            % Hold on the trace of the modal participation
    ylim([-autorange*1.5 autorange*1.5]) 
    if i == 1
        xlabel('Mode Index')
        ylabel('Amplitude')
        title('Snapshot of the Modal Participation')
        pbaspect([1.618 1 1])
        xticks([1 2 3 4 5 6 7 8 9 10 11])
        xticklabels({'Base','#1','#2','#3','#4','#5','#6','#7','#8','#9', '#10'})
    end
end
disp('Modal Participation Visualization Done!')
stem([0 q_am(i,1:n)],'r')                              % Retain the last played position

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
