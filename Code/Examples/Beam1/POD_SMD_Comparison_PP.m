%% SMOOTH MODE DECOMPOSITION
% ------------------------------------------------------------------------
% SMD of a spatial-temporally underdetermined problem cannot resolve unique
% solution of the system
% ------------------------------------------------------------------------
% Set up the environment
clear
% Environment Setup
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
% Load the data
load('Data.mat');
Y_raw = Data.Y;
dx = Data.SpatialResolution;
fs = Data.fs;
Phi_x = Data.Phi_x;
Phi_am = real(Phi_x);
[m, n] = size(Y_raw); % The size of the raw data
Eta = real(Data.Eta(:,1:2:end));
%% Smooth Coordinate Decomposition
PrintFlag = 1;
Min_rr_SMD = m/n;
% Set up the resampling rate 
rr = 1; % TO GENERATE THE ONE USED IN THE REPORT rr = 460;
var_Y_mean = mean(var(Y_raw));
SNR = inf;
NoiseLevel = var_Y_mean/SNR;
% if rr < Min_rr_SMD
%     error(['The minimum size of resampling rate has to be larger than:',num2str(Min_rr_SMD)])
% end
rng(1)
% Data Resampling (down-sampling)
Y_raw_noise = Y_raw + NoiseLevel*randn(size(Y_raw));
Y = Y_raw_noise(1:rr:end, :); 

figure(1),clf
imagesc(Y)
view(270,-90)
pbaspect([1 2.5 1])
xlabel('Time Samples')
ylabel('Spatial Samples')
set(gcf,'papersize', [6 2.5])
set(gcf,'paperposition', [0 0 6 2.5])
mycolormap = mycolorbar('Viridis');
colormap(mycolormap)

DY = fs*diff(Y);
[soc, sov, spm, som, S1, S2, U1, U2] = sod(Y, DY, ' ', 'true');
[poc, pov, pom] = svd(Y, 'econ');

% Apply the modal assurance criterion to obtain the sorting indices
[~, MIndx_scd, SIndx_scd] = MAC(som, real(Phi_x));
[~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));

% Plot estimation results
figure(2),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    pom_n = normalize(pom(:,SIndx_pod(i)),'norm');
    som_n = normalize(som(:,SIndx_scd(i)),'norm');
    Phi_am_n = normalize(Phi_am(:,i),'norm');
    % Check the orientation of the modes
    dir_pod = sign(pom_n'*Phi_am_n);
    % Check the orientation of the modes
    dir_sod = sign(som_n'*Phi_am_n);
    plot([0; dir_pod*pom_n])
    hold on

    % SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS
    plot([0;dir_sod*som_n])
    plot([0; Phi_am_n], 'k--')
    axis tight

    ylabel(['Mode ',num2str(i)])

    if count == 1
        legend(['$\widehat\phi_\mathrm{pod}$'],['$\widehat\phi_\mathrm{smd}$'],['$\widehat\phi$'],'NumColumns',3, 'location','northoutside')
    end

    if count == 9 || count == 10
        xlabel('Noal Point')
    end
    count = count + 1;
    ylim([-0.15 0.15])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])

sgtitle('Mode Shape Estiamtion between POD and SCD')

% 
% figure(2),clf
% plot(normalize(soc(:,1),'norm'))
% hold on
% plot(normalize(Eta(:,1),'norm'))

figure(3),clf,MAC(pom, real(Phi_x));
sgtitle('Modal Assurance Criteria between $\Phi_{POD}$ and $\Phi$')
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

figure(4),clf
MAC(som, real(Phi_x));
%     MAC(som_smd, real(Phi_x(:,1:size(som_smd,2))));
sgtitle(['Modal Assurance Criteria between $\Phi_{SCD}$ and $\Phi$; rr = ',num2str(rr)])
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

if PrintFlag == 1
    figure(1), print(['Y_resampled_rr_',num2str(rr),'.png'],'-dpng','-r600')
    figure(2), print('POD_SCD_comparison_pp_beam.png','-dpng','-r600')
    figure(3), print(['MAC_POD_rr_',num2str(rr),'.png'],'-dpng','-r600')
    figure(4), print(['MAC_SCD_rr_',num2str(rr),'.png'],'-dpng','-r600')
end

%=============== Temporally underdetermined problem =====================
%% Temporal Dimensionality Reduction (TDR-1) --- TSVD based
close all

% 1. Direct POD
[poc, pov, pom] = svd(Y_raw, 'econ');
plot(diag(pov),'--*')
set(gca,'yscale','log')
pbaspect([2.5 1 1])
xlabel('Index')
ylabel('POV')

% 2. TSVD based temporal dimensionality reduction
r = 20;
Q = poc(:, 1:r);
Yr = Q'*Y;
delYr = GenFiniteDiff(Yr', dx, 'c2')'; % Center finite difference

% 3. TSMD
[som_tsmd, sov_tsmd, spm_tsmd, sc_tsmd, S1_tsmd, S2_tsmd, U_tsmd, V_tsmd] = sod(Yr', delYr');

[m_tsmd, n_tsmd] = size(som_tsmd);

% 4. Transform back to the origina higher dim.
sc = Q*sc_tsmd;

% Apply the modal assurance criterion to obtain the sorting indices
[~, MIndx_tsmd, SIndx_tsmd] = MAC(som_tsmd, real(Phi_x));
[~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));
% Apply the modal assurance criterion to obtain the sorting indices
% [~, MIndx_tsmd, SIndx_tsmd] = MAC(sc, Eta);
% [~, MIndx_pod, SIndx_pod] = MAC(poc, Eta);

Phi_am = real(Phi_x);

% Plot estimation results
figure(2),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    pom_n = normalize(pom(:,SIndx_pod(i)),'norm');
    som_n = normalize(som_tsmd(:,SIndx_tsmd(i)),'norm');
    Phi_am_n = normalize(Phi_am(:,i),'norm');
    % Check the orientation of the modes
    dir_pod = sign(pom_n'*Phi_am_n);
    % Check the orientation of the modes
    dir_sod = sign(som_n'*Phi_am_n);
    plot([0; dir_pod*pom_n])
    hold on

    % SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS
    plot([0;dir_sod*som_n])
    plot([0; Phi_am_n], 'k--')
    axis tight

    ylabel(['Mode ',num2str(i)])

    if count == 1
        legend(['$\widehat\phi_\mathrm{pod}$'],['$\widehat\phi_\mathrm{smd}$'],['$\phi$'],'NumColumns',3, 'location','northoutside')
    end

    if count == 9 || count == 10
        xlabel('Noal Point')
    end
    count = count + 1;
    ylim([-0.15 0.15])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])

sgt = sgtitle(['Mode Shape Estiamtion between POD and TSMD-P; $rr = ',num2str(rr),';r = ',num2str(r),'$']);
sgt.FontSize = 12;

figure(3),clf,MAC(pom, real(Phi_x));
sgtitle('Modal Assurance Criteria between $\Phi_{POD}$ and $\Phi$')
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

figure(4),clf
MAC(som_tsmd, real(Phi_x));
%     MAC(som_smd, real(Phi_x(:,1:size(som_smd,2))));
sgt = sgtitle(['Modal Assurance Criteria between $\Phi_{TSMD-P}$ and $\Phi; rr = ',num2str(rr),';r = ',num2str(r),'$']);
sgt.FontSize = 12;

subplot(121)
xlim([0.5 r + 0.5])
ylim([0.5 r + 0.5])
subplot(122)
xlim([0.5 r + 0.5])
ylim([0.5 r + 0.5])

set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

% Plot estimation results
figure(5),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    poc_n = normalize(poc(:,SIndx_pod(i)),'norm');
    sc_n = normalize(sc(:,SIndx_tsmd(i)),'norm');
    Eta_n = normalize(Eta(:,i),'norm');
    % Check the orientation of the modes
    dir_poc = sign(poc_n'*Eta_n);
    % Check the orientation of the modes
    dir_sc = sign(sc_n'*Eta_n);

    [pxx_poc, fxx] = pwelch(dir_poc*poc_n, boxcar(m/4), m/8, [], fs);
    [pxx_sc, ~] = pwelch(dir_sc*sc_n, boxcar(m/4), m/8, [], fs);
    [pxx_eta, ~] = pwelch(Eta_n, boxcar(m/4), m/8, [], fs);

        plot([dir_poc*poc_n])
%     plot(fxx, 10*log10(pxx_poc))
    hold on

    % SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS
        plot(dir_sc*sc_n)
%     plot(fxx, 10*log10(pxx_sc))
        plot(Eta_n,'--k')
%     plot(fxx, 10*log10(pxx_eta), 'k--')
    axis tight

    ylabel(['Mode ',num2str(i)])

    if count == 1
        legend(['$\widehat\eta_\mathrm{pod}$'],['$\widehat\eta_\mathrm{tsmd}$'],['$\eta$'],'NumColumns',3, 'location','northoutside')
    end

    if count == 9 || count == 10
        xlabel('Sample Time')
    end
    count = count + 1;
%     xlim([0 5])
%     ylim([-0.02 0.02])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])
sgt = sgtitle(['Modal Coordinate Estiamtion between POD and TSMD-P; $rr = ',num2str(rr),';r = ',num2str(r),'$']);
sgt.FontSize = 12;

% 
if PrintFlag == 1
    figure(2), print('POD_TSMDP_phi_comparison_pp_beam.png','-dpng','-r600')
%     figure(3), print(['MAC_POD_rr_',num2str(rr),'.png'],'-dpng','-r600')
    figure(4), print(['MAC_TSMDP_rr_',num2str(rr),'_r_',num2str(r),'.png'],'-dpng','-r600')
    figure(5), print('POD_TSMDP_eta_comparison_pp_beam.png','-dpng','-r600')
end
% figure(4)
% print('POD_SMD_comparison_pp_beam.png','-dpng','-r600')
% print('full_POD_SMD_r100_comparison_pp_beam.png','-dpng','-r600')
%% Temporal Dimensionality Reduction (TDR-2) --- TSVD-D based
PrintFlag = 1;
Min_rr_SMD = m/n;
% Set up the resampling rate 
rr = 1; % TO GENERATE THE ONE USED IN THE REPORT rr = 460;
var_Y_mean = mean(var(Y_raw));
SNR = inf;
NoiseLevel = var_Y_mean/SNR;
% if rr < Min_rr_SMD
%     error(['The minimum size of resampling rate has to be larger than:',num2str(Min_rr_SMD)])
% end
rng(1)
% Data Resampling (down-sampling)
Y_raw_noise = Y_raw + NoiseLevel*randn(size(Y_raw));
Y = Y_raw_noise(1:rr:end, :); 

% 1. Direct POD to concatenated matrices [Y, DY]
[poc, pov, pom] = svd(Y_raw, 'econ');
DY = fs*diff(Y);
Y = Y(1:end-1,:);
[Ut, St, Vt] = svd([Y,DY], 'econ');
delY = GenFiniteDiff(Y', dx, 'c2')'; % Center finite difference

% 2. TSVD based temporal dimensionality reduction
r = 20;
Q = Ut(:, 1:r);
Yr = Q'*Y;
delYr = Q'*delY; % Center finite difference

% 3. TSMD
[som_tsmd, sov_tsmd, spm_tsmd, sc_tsmd, S1_tsmd, S2_tsmd, U_tsmd, V_tsmd] = sod(Yr', delYr');

[m_tsmd, n_tsmd] = size(som_tsmd);

% 4. Transform back to the origina higher dim.
sc = Q*sc_tsmd;

% Apply the modal assurance criterion to obtain the sorting indices
[~, MIndx_tsmd, SIndx_tsmd] = MAC(som_tsmd, real(Phi_x));
[~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));

Phi_am = real(Phi_x);

% Plot estimation results
figure(2),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    pom_n = normalize(pom(:,SIndx_pod(i)),'norm');
    som_n = normalize(som_tsmd(:,SIndx_tsmd(i)),'norm');
    Phi_am_n = normalize(Phi_am(:,i),'norm');
    % Check the orientation of the modes
    dir_pod = sign(pom_n'*Phi_am_n);
    % Check the orientation of the modes
    dir_sod = sign(som_n'*Phi_am_n);
    plot([0; dir_pod*pom_n])
    hold on

    % SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS
    plot([0;dir_sod*som_n])
    plot([0; Phi_am_n], 'k--')
    axis tight

    ylabel(['Mode ',num2str(i)])

    if count == 1
        legend(['$\widehat\phi_\mathrm{pod}$'],['$\widehat\phi_\mathrm{smd}$'],['$\phi$'],'NumColumns',3, 'location','northoutside')
    end

    if count == 9 || count == 10
        xlabel('Noal Point')
    end
    count = count + 1;
    ylim([-0.15 0.15])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])

sgt = sgtitle(['Mode Shape Estiamtion between POD and TSMD-D; $rr = ',num2str(rr),';r = ',num2str(r),'$']);
sgt.FontSize = 12;

figure(3),clf,MAC(pom, real(Phi_x));
sgtitle('Modal Assurance Criteria between $\Phi_{POD}$ and $\Phi$')
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

figure(4),clf
MAC(som_tsmd, real(Phi_x));
%     MAC(som_smd, real(Phi_x(:,1:size(som_smd,2))));
sgt = sgtitle(['Modal Assurance Criteria between $\Phi_{TSMD-D}$ and $\Phi; rr = ',num2str(rr),';r = ',num2str(r),'$']);
sgt.FontSize = 12;

subplot(121)
xlim([0.5 r + 0.5])
ylim([0.5 r + 0.5])
subplot(122)
xlim([0.5 r + 0.5])
ylim([0.5 r + 0.5])

set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

% Plot estimation results
figure(5),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    poc_n = normalize(poc(1:end-1,SIndx_pod(i)),'norm');
    sc_n = normalize(sc(:,SIndx_tsmd(i)),'norm');
    Eta_n = normalize(Eta(1:end-1,i),'norm');
    % Check the orientation of the modes
    dir_poc = sign(poc_n'*Eta_n);
    % Check the orientation of the modes
    dir_sc = sign(sc_n'*Eta_n);

    [pxx_poc, fxx] = pwelch(dir_poc*poc_n, boxcar(m/4), m/8, [], fs);
    [pxx_sc, ~] = pwelch(dir_sc*sc_n, boxcar(m/4), m/8, [], fs);
    [pxx_eta, ~] = pwelch(Eta_n, boxcar(m/4), m/8, [], fs);

        plot([dir_poc*poc_n])
%     plot(fxx, 10*log10(pxx_poc))
    hold on

    % SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS
        plot(dir_sc*sc_n)
%     plot(fxx, 10*log10(pxx_sc))
        plot(Eta_n,'--k')
%     plot(fxx, 10*log10(pxx_eta), 'k--')
    axis tight

    ylabel(['Mode ',num2str(i)])

    if count == 1
        legend(['$\widehat\eta_\mathrm{pod}$'],['$\widehat\eta_\mathrm{tsmd}$'],['$\eta$'],'NumColumns',3, 'location','northoutside')
    end

    if count == 9 || count == 10
        xlabel('Sample Time')
    end
    count = count + 1;
%     xlim([0 5])
%     ylim([-0.02 0.02])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])
sgt = sgtitle(['Modal Coordinate Estiamtion between POD and TSMD-D; $rr = ',num2str(rr),';r = ',num2str(r),'$']);
sgt.FontSize = 12;

% 
% if PrintFlag == 1
    figure(2), print('POD_TSMDD_phi_comparison_pp_beam.png','-dpng','-r600')
% %     figure(3), print(['MAC_POD_rr_',num2str(rr),'.png'],'-dpng','-r600')
    figure(4), print(['MAC_TSMDD_rr_',num2str(rr),'_r_',num2str(r),'.png'],'-dpng','-r600')
    figure(5), print('POD_TSMDD_eta_comparison_pp_beam.png','-dpng','-r600')
% end

%% Temporal Dimensionality Reduction (TDR-3) --- TSVD-DEL based

% 1. Direct POD to concatenated matrices [Y, DY]
[poc, pov, pom] = svd(Y_raw, 'econ');

% DY = fs*diff(Y);
delY = GenFiniteDiff(Y', dx, 'c2')'; % Center finite difference
[Ut, St, Vt] = svd([Y,delY], 'econ');


% 2. TSVD based temporal dimensionality reduction
r = 20;
Q = Ut(:, 1:r);
Yr = Q'*Y;
delYr = Q'*delY; % Center finite difference

% 3. TSMD
[som_tsmd, sov_tsmd, spm_tsmd, sc_tsmd, S1_tsmd, S2_tsmd, U_tsmd, V_tsmd] = sod(Yr', delYr');

[m_tsmd, n_tsmd] = size(som_tsmd);

% 4. Transform back to the origina higher dim.
sc = Q*sc_tsmd;

% Apply the modal assurance criterion to obtain the sorting indices
[~, MIndx_tsmd, SIndx_tsmd] = MAC(som_tsmd, real(Phi_x));
[~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));

Phi_am = real(Phi_x);

% Plot estimation results
figure(2),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    pom_n = normalize(pom(:,SIndx_pod(i)),'norm');
    som_n = normalize(som_tsmd(:,SIndx_tsmd(i)),'norm');
    Phi_am_n = normalize(Phi_am(:,i),'norm');
    % Check the orientation of the modes
    dir_pod = sign(pom_n'*Phi_am_n);
    % Check the orientation of the modes
    dir_sod = sign(som_n'*Phi_am_n);
    plot([0; dir_pod*pom_n])
    hold on

    % SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS
    plot([0;dir_sod*som_n])
    plot([0; Phi_am_n], 'k--')
    axis tight

    ylabel(['Mode ',num2str(i)])

    if count == 1
        legend(['$\widehat\phi_\mathrm{pod}$'],['$\widehat\phi_\mathrm{smd}$'],['$\phi$'],'NumColumns',3, 'location','northoutside')
    end

    if count == 9 || count == 10
        xlabel('Noal Point')
    end
    count = count + 1;
    ylim([-0.15 0.15])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])

sgt = sgtitle(['Mode Shape Estiamtion between POD and TSMD-N; $rr = ',num2str(rr),';r = ',num2str(r),'$']);
sgt.FontSize = 12;

figure(3),clf,MAC(pom, real(Phi_x));
sgtitle('Modal Assurance Criteria between $\Phi_{POD}$ and $\Phi$')
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

figure(4),clf
MAC(som_tsmd, real(Phi_x));
%     MAC(som_smd, real(Phi_x(:,1:size(som_smd,2))));
sgt = sgtitle(['Modal Assurance Criteria between $\Phi_{TSMD-N}$ and $\Phi; rr = ',num2str(rr),';r = ',num2str(r),'$']);
sgt.FontSize = 12;

subplot(121)
xlim([0.5 r + 0.5])
ylim([0.5 r + 0.5])
subplot(122)
xlim([0.5 r + 0.5])
ylim([0.5 r + 0.5])

set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

% Plot estimation results
figure(5),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    poc_n = normalize(poc(1:end-1,SIndx_pod(i)),'norm');
    sc_n = normalize(sc(:,SIndx_tsmd(i)),'norm');
    Eta_n = normalize(Eta(1:end-1,i),'norm');
    % Check the orientation of the modes
    dir_poc = sign(poc_n'*Eta_n);
    % Check the orientation of the modes
    dir_sc = sign(sc_n'*Eta_n);

    [pxx_poc, fxx] = pwelch(dir_poc*poc_n, boxcar(m/4), m/8, [], fs);
    [pxx_sc, ~] = pwelch(dir_sc*sc_n, boxcar(m/4), m/8, [], fs);
    [pxx_eta, ~] = pwelch(Eta_n, boxcar(m/4), m/8, [], fs);

        plot([dir_poc*poc_n])
%     plot(fxx, 10*log10(pxx_poc))
    hold on

    % SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS
        plot(dir_sc*sc_n)
%     plot(fxx, 10*log10(pxx_sc))
        plot(Eta_n,'--k')
%     plot(fxx, 10*log10(pxx_eta), 'k--')
    axis tight

    ylabel(['Mode ',num2str(i)])

    if count == 1
        legend(['$\widehat\eta_\mathrm{pod}$'],['$\widehat\eta_\mathrm{tsmd}$'],['$\eta$'],'NumColumns',3, 'location','northoutside')
    end

    if count == 9 || count == 10
        xlabel('Sample Time')
    end
    count = count + 1;
%     xlim([0 5])
%     ylim([-0.02 0.02])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])
sgt = sgtitle(['Modal Coordinate Estiamtion between POD and TSMD-N; $rr = ',num2str(rr),';r = ',num2str(r),'$']);
sgt.FontSize = 12;


if PrintFlag == 1
    figure(2), print('POD_TSMDN_phi_comparison_pp_beam.png','-dpng','-r600')
    figure(4), print(['MAC_TSMDN_rr_',num2str(rr),'_r_',num2str(r),'.png'],'-dpng','-r600')
    figure(5), print('POD_TSMDN_eta_comparison_pp_beam.png','-dpng','-r600')
end

%% Grid search for the resamping cretiria
% ------------------------------------------------------------------------
%           INVESTIGATE HOW RR AFFECT THE SMD ACCURACY
% ------------------------------------------------------------------------
N = 600;                     % Total number of grids
% Allocate memory
error_smd = zeros(N,1);      % Error from smd 
error_pod = zeros(N,1);      % Error from pod
% Set up environment
addpath('L:\My Drive\Graduate study\Research\Projects\OS')
% Grid search begins
for j = 1:N
    progress_bar(j, N, 'Grid Searching Minimum Error Resampling Rate')
    % SMD
    addpath('L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox')
    % addpath('G:\LNN_NNM_ID\Code')
    Y = Y_raw(1:50+j:end, :); % 188， 210， 267, 310, 400, 464, 479
    % Y = Y + 0.0001*randn(size(Y));
    delY = GenFiniteDiff(Y', dx, 'c2')';

    % DIRECT SMD
    [som_smd, sov_smd, spm_smd, soc_smd, S1_smd, S2_smd, U_smd, V_smd] = sod(Y', delY');

    [poc, pov, pom] = svd(Y, 'econ');

    % Shrink the size in the spatial domain
    [C, MIndx_smd, SIndx_smd] = MAC(som_smd, real(Phi_x));
    [~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));
    % Plot estimation results
    count = 1;
    strt = 1;
    Nmodes = 20;   % Number of modes to consider in the cumulative error calc
    error_smd_temp = zeros(Nmodes,1);
    error_pod_temp = zeros(Nmodes,1);

    for i = 1:Nmodes
        % Get the normalized modes (according to some normalization scheme)
        pom_n = normalize(pom(:,SIndx_pod(i)),'norm');
        som_n = normalize(som_smd(:,SIndx_smd(i)),'norm');
        Phi_am_n = normalize(Phi_am(:,i),'norm');
        % Check the orientation of the modes
        dir_pod = sign(pom_n'*Phi_am_n);
        % Check the orientation of the modes
        dir_sod = sign(som_n'*Phi_am_n);
        hold on
        count = count + 1;
        % COMPUTE THE CUMULATIVE SQUARE ERROR
        error_smd_temp(i) = sum((Phi_am_n - dir_sod*som_n).^2);
        error_pod_temp(i) = sum((Phi_am_n - dir_pod*pom_n).^2);
    end
    % SUM OVER ONE CASE 
    error_smd(j) = sum(error_smd_temp);
    error_pod(j) = sum(error_pod_temp);
end
% Find mininum estimated cumulative error
[min_error_smd, Indx_smd] = min(error_smd);
[min_error_pod, Indx_pod] = min(error_pod);
% Plot the results
figure(7),clf
plot(error_smd)
hold on
plot(error_pod)
plot(Indx_smd, min_error_smd,'ro')
plot(Indx_pod, min_error_pod,'rx')
xticks([0:50:600])
xticklabels({'50','100','150','200','250','300','350','400','450','500','550','600','650'})
xlabel('Resampling Rate - $r$')
ylabel('Cumulative Error - $\sum_{i = 1}^{20}(\hat\phi_i - \phi_i)^2$')
grid on
set(gcf,'papersize',[6 2.5])
set(gcf,'paperposition',[0 0 6 2.5])
pbaspect([2.5 1 1])
legend('SMD','POD')
print('CumError_POD_SMD_r.png','-dpng','-r600')
%% TSMD
% SMD
[poc, pov, pom] = svd(Y_raw,'econ');
[Uc, Sc, Vc] = svd([Y; delY]);
plot(diag(Sc))
r = 120; %110
Vr = Vc(:,1:r);  % New basis for the data
Yt = Y*Vr;      % Project the data onto the new basis
delYt = delY*Vr;    % Project the derivative onto the same basis

addpath('L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox')
tic
[somt, sovt, spmt, soct, S1t, S2t, Ut, Vt] = sod(Yt', delYt');
som_rsmd = Vr*somt;
toc

% Shrink the size in the spatial domain
[C, MIndx_rsmd, SIndx_rsmd] = MAC(som_rsmd, real(Phi_x));
[~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));

figure(5), clf
count = 1;
for i = 1:10
subplot(5,2,count)
% Get the normalized modes (according to some normalization scheme)
pom_n = normalize(pom(:,SIndx_pod(i)),'norm');
som_n = normalize(som_rsmd(:,SIndx_rsmd(i)),'norm');
Phi_am_n = normalize(real(Phi_x(:,i)),'norm');
% Check the orientation of the modes
dir_pod = sign(pom_n'*Phi_am_n);
% Check the orientation of the modes
dir_sod = sign(som_n'*Phi_am_n);
plot([0; dir_pod*pom_n])
hold on

% SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS
plot([0;dir_sod*som_n])
plot([0; Phi_am_n], 'k--')
axis tight
ylabel([-0.07 0.07])
ylabel(['Mode ',num2str(i)])

if count == 1
legend(['$\widehat\phi_\mathrm{pod}$'],['$\widehat\phi_\mathrm{rsmd}$'],['$\widehat\phi$'],'NumColumns',3, 'location','northoutside')
end

if count == 9 || count == 10
    xlabel('Noal Point')
end
count = count + 1;
end
sgtitle()
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])
sgtitle('Mode Shape Estiamtion between POD and SOD (RSMD)')
print(['POD_RSMD_comparison_pp_beam_r_',num2str(r),'.png'],'-dpng','-r600')

% sgtitle('Mode Shape Estiamtion between POD and SOD (RSMD)')
% Shrink the size in the spatial domain
figure(6),clf
MAC(som_rsmd, real(Phi_x))
sgtitle('Modal Assurance Criteria between $\Phi_{RSMD}$ and $\Phi$')
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])
print(['MAC_POD_RSMD_r_',num2str(r),'.png'],'-dpng','-r600')


figure(7),clf
MAC(pom, real(Phi_x))
sgtitle('Modal Assurance Criteria between $\Phi_{POD}$ and $\Phi$')
%% Grid search for the resamping cretiria
% ------------------------------------------------------------------------
%           INVESTIGATE HOW RR AFFECT THE RSMD ACCURACY
%          IN COMPARISON TO THE FULL POD-RESOLVED MODES
% ------------------------------------------------------------------------
r = size(Y,1)+1:1:round(size(Y,2)); 
% SVD to the resampled data
[Uc, Sc, Vc] = svd([Y; delY]);
[~, N] = size(r);
% Allocate memory
error_rsmd = zeros(N,1);      % Error from smd 
error_pod = zeros(N,1);      % Error from pod
% Set up environment
addpath('L:\My Drive\Graduate study\Research\Projects\OS')
% Grid search begins
for j = 1:N
    progress_bar(j, N, 'Grid Searching Minimum Error Resampling Rate')
    % SMD
    addpath('L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox')
    
    plot(diag(Sc))

    Vr = Vc(:,1:r(j));  % New basis for the data
    Yt = Y*Vr;      % Project the data onto the new basis
    delYt = delY*Vr;    % Project the derivative onto the same basis
   
    % Smooth mode decomposition
    [somt, sovt, spmt, soct, S1t, S2t, Ut, Vt] = sod(Yt', delYt');
    som_rsmd = Vr*somt;
    
    % Shrink the size in the spatial domain
    [C, MIndx_rsmd, SIndx_rsmd] = MAC(som_rsmd, real(Phi_x));
    [~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));
    
    % Plot estimation results
    count = 1;
    strt = 1;
    Nmodes = 10;   % Number of modes to consider in the cumulative error calc
    error_rsmd_temp = zeros(Nmodes,1);
    error_pod_temp = zeros(Nmodes,1);

    for i = 1:Nmodes
        % Get the normalized modes (according to some normalization scheme)
        pom_n = normalize(pom(:,SIndx_pod(i)),'norm');
        som_n = normalize(som_rsmd(:,SIndx_rsmd(i)),'norm');
        Phi_am_n = normalize(Phi_am(:,i),'norm');
        % Check the orientation of the modes
        dir_pod = sign(pom_n'*Phi_am_n);
        % Check the orientation of the modes
        dir_sod = sign(som_n'*Phi_am_n);
        hold on
        count = count + 1;
        % COMPUTE THE CUMULATIVE SQUARE ERROR
        error_rsmd_temp(i) = sum((Phi_am_n - dir_sod*som_n).^2);
        error_pod_temp(i) = sum((Phi_am_n - dir_pod*pom_n).^2);
    end
    % SUM OVER ONE CASE 
    error_rsmd(j) = sum(error_rsmd_temp);
    error_pod(j) = sum(error_pod_temp);
end
% Find mininum estimated cumulative error
[min_error_rsmd, Indx_rsmd] = min(error_rsmd);
[min_error_pod, Indx_pod] = min(error_pod);
% Plot the results
figure(7),clf
plot(error_rsmd)
hold on
plot(error_pod)
plot(Indx_rsmd, min_error_rsmd,'ro')
plot(Indx_pod, min_error_pod,'rx')
xticks([1:50:r(end-1)])

xticklabels(cellstr(string(r(floor(1:50:length(r))))))
xlabel('Size of the reduced spatial space - $r$')
ylabel(['Cumulative Error - $\sum_{i = 1}^{',num2str(Nmodes),'}(\hat\phi_i - \phi_i)^2$'])
grid on
set(gcf,'papersize',[6 2.5])
set(gcf,'paperposition',[0 0 6 2.5])
pbaspect([2.5 1 1])
legend('RSMD','POD')
axis tight
% print('CumError_POD_RSMD_r.png','-dpng','-r600')
%% TSMD
Min_rr_SMD = m/n;
% Set up the resampling rate 
rr = 1; % TO GENERATE THE ONE USED IN THE REPORT rr = 460;
var_Y_mean = mean(var(Y_raw));
SNR = inf;
NoiseLevel = var_Y_mean/SNR;
rng(1)
% Data Resampling (down-sampling)
Y_raw_noise = Y_raw + NoiseLevel*randn(size(Y_raw));
Y = Y_raw_noise(1:rr:end, :); 

% Apply finite difference method to the data
delY = GenFiniteDiff(Y', dx, 'c2')'; % Center difference
dY = diff(Y);
% SMD
[poc, pov, pom] = svd(Y,'econ');
[Uc, Sc, Vc] = svd([Y, delY], 'econ');

plot(diag(Sc))
r = 10;
Ur = Uc(:,1:r);  % New basis for the data
Yt = Ur'*Y;      % Project the data onto the new basis
delYt = Ur'*delY;    % Project the derivative onto the same basis
addpath('L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox')

% Conduct SMD
[somt, sovt, spmt, soct, S1t, S2t, Ut, Vt] = sod(Yt', delYt');
soc_tsmd = Ur*soct;

% Shrink the size in the spatial domain
[C, MIndx_rsmd, SIndx_rsmd] = MAC(som_rsmd, real(Phi_x));
[~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));

%%
r = 100;
pocr = poc(:,1:r);
Yt = pocr'*Y;
delYt = GenFiniteDiff(Yt', dx, 'c2')'; % Center difference
addpath('L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox')

% Conduct SMD
[somt, sovt, spmt, soct, S1t, S2t, Ut, Vt] = sod(Yt', delYt');
soc_tsmd = pocr*soct;

plot(soc_tsmd(:,1:10))

% Shrink the size in the spatial domain
[C, MIndx_rsmd, SIndx_rsmd] = MAC(som_tsmd, real(Phi_x));
[~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));

%% TSMD
Min_rr_SMD = m/n;
% Set up the resampling rate 
rr = 1; % TO GENERATE THE ONE USED IN THE REPORT rr = 460;
var_Y_mean = mean(var(Y_raw));
SNR = inf;
NoiseLevel = var_Y_mean/SNR;
rng(1)
% Data Resampling (down-sampling)
Y_raw_noise = Y_raw + NoiseLevel*randn(size(Y_raw));
Y = Y_raw_noise(1:rr:end, :); 

% Apply finite difference method to the data
dY = diff(Y);
% SMD
[poc, pov, pom] = svd(Y,'econ');
[Uc, Sc, Vc] = svd([Y(1:end-1,:), dY], 'econ');

plot(diag(Sc))
r = 20;
Ur = Uc(:,1:r);             % New basis for the data
Yt = Ur'*Y(1:end-1,:);      % Project the data onto the new basis
delYt = Ur'*dY;             % Project the derivative onto the same basis
addpath('L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox')

% Conduct SMD
[somt, sovt, spmt, soct, S1t, S2t, Ut, Vt] = sod(Yt', delYt');
soc_tsmd = Ur*soct;

plot(soc_tsmd(:,3))
plot(Ut(:,1:5))
% Shrink the size in the spatial domain
[C, MIndx_rsmd, SIndx_rsmd] = MAC(som_rsmd, real(Phi_x));
[~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));