%% Case IV - N10n10201m601 - 1
% ------------------------------------------------------------------------
% SMD of a spatial-temporally underdetermined problem cannot resolve unique
% solution of the system
% ------------------------------------------------------------------------
% Set up the environment
clear
close all
addpath('L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox')
addpath('L:\My Drive\Graduate study\Research\Projects\OS')
set(0,'DefaultFigureWindowStyle','docked')
[Markers, LineStyles, Mycolors] = mystyle();
DefColors = Mycolors.defaultcolor;
load('Data.mat');

CurrentCase = split(cd,'\');
CurrentCase = CurrentCase{end};
%% Data Assignment and truncation
mnew = 1200;
rr = 1;
Y_raw = Data.Y;
TY_raw = Data.TY;
a0 = Data.Dimensions.a0;
b0 = Data.Dimensions.b0;
% matrix is full row rank
Y_raw = Y_raw(1:rr:mnew,:);
TY_raw = TY_raw(:,:,1:rr:mnew);
dx = Data.SpatialResolution.dx;
dy = Data.SpatialResolution.dy;
fs = Data.fs;
dt = 1/fs;
tend = Data.tend;
t = 0:dt:tend;
Phi_x = Data.Phi_x;
Phi_am = real(Phi_x);
[m, n] = size(Y_raw); % The size of the raw data
[mm, nn, ~] = size(TY_raw);
Eta = real(Data.Eta(1:rr:mnew,1:2:end));
t = t(1:rr:mnew);
%% Visualize the True Modes and Coordinates
windsize = mnew/2;
noverlap = windsize/8;
xwindsize = size(Phi_x,1)/2;
xnoverlap = xwindsize/4;
nrows = 10;
figure(100),clf
for i = 1:nrows
    subplot(nrows,4,4*(i-1)+1)
    imagesc(reshape(normalize(real(Phi_x(:,i)),'norm'),mm,nn))
    mycolorbar()
    if i == nrows
        xlabel('x')
    end
    ylabel('y')
    if i == 1
        title('Modes --- $\Phi_i$')
    end
    subplot(nrows,4,4*(i-1)+2)
    plot(t,Eta(:,i))
    if i == nrows
        xlabel('Time (seconds)')
    end
    ylabel(['$q_',num2str(i),'(t)$'])
    axis tight
    if i == 1
        title('Coords. --- $q_i$')
    end

    subplot(nrows,4,4*(i-1)+3)
    [peta, fxx] = pwelch(Eta,hanning(windsize),noverlap,[],fs);
    plot(fxx, 10*log10(peta),'color','#B6B6B4')
    hold on
    plot(fxx, 10*log10(peta(:,i)),'color', DefColors{1})
    if i == nrows
        xlabel('Frequency (Hz)')
    end
    ylabel('PSD (Power/Hz)')
    xlim([0 4])
    if i == 1
        title('PSD[$q_i(t)$]')
    end

    subplot(nrows,4,4*(i-1)+4)
    [px, fxx] = pwelch(real(Phi_x),hanning(xwindsize),xnoverlap);
    plot(fxx, 10*log10(px),'color','#B6B6B4')
    hold on
    plot(fxx, 10*log10(px(:,i)),'color', DefColors{2})
%     imagesc(10*log10(abs(fftshift(fft2(reshape(real(Phi_x(:,1)),mm,nn)))).^2));

    if i == nrows
        xlabel('Frequency (Hz)')
    end
    ylabel('PSD (Power/Hz)')
    if i == 1
        title('PSD[$\phi_i$]')
    end
end

set(gcf,'papersize',[10 14])
set(gcf,'paperposition',[0 0 10 14])
print(['Phi_info_',CurrentCase,'.png'],'-dpng','-r600')
%% Smooth Mode Decomposition
PrintFlag = 0;
SpectralPlot = 1;
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
% Put data into a tensor
TY = Y2TY(Y, mm, nn);
[delY, ~] = TY2delY(TY, dx, dy);

% IT WAS A WRONG IMPLIMENTATION HERE!!!!
% delY = GenFiniteDiff(Y','c2');

% ------------------- VISUALIZATION OF RESAMPLED DATA --------------------
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
% Direct POD
[poc, pov, pom] = svd(Y, 'econ');

plot((diag(pov)))
set(gca,'yscale','log')
xlabel('Index of POM')
ylabel('POV')
pbaspect([2 1 1])
% DIRECT SMD
tic
[som_smd, sov_smd, spm_smd, soc_smd, S1_smd, S2_smd, U_smd, V_smd] = sod(Y', delY',' ','False');
toc
% Apply the modal assurance criterion to obtain the sorting indices
[~, MIndx_smd, SIndx_smd] = MAC(som_smd, real(Phi_x));
[~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));


% ------------------- VISUALIZATION OF ESTIMATED PHI --------------------
% Plot estimation results
figure(2),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    pom_n = normalize(pom(:,SIndx_pod(i)),'norm');
    som_n = normalize(som_smd(:,SIndx_smd(i)),'norm');
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
    ylim([-0.05 0.05])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])

sgtitle('Mode Shape Estiamtion between POD and SCD')



% Plot the comparison between POD and SMD
nrows = 8;
strt = 1;
figure(3),clf
for i = 1:nrows
subplot(nrows,2,2*(i-1)+1)
imagesc(vec2snap(som_smd(:,SIndx_smd(i)),mm,nn))
pbaspect([b0 a0 1])
subplot(nrows,2,2*(i-1)+2)
plot(soc_smd(:,SIndx_smd(i)))
pbaspect([b0 a0 1])
end

nrows = 8;
strt = 1;
figure(4),clf
for i = 1:nrows
subplot(nrows,2,2*(i-1)+1)
imagesc(vec2snap(pom(:,strt+i-1),mm,nn))
pbaspect([b0 a0 1])
subplot(nrows,2,2*(i-1)+2)
plot(poc(:,strt+i-1))
pbaspect([b0 a0 1])
end

figure(5),clf,MAC(pom, real(Phi_x));
sgtitle('Modal Assurance Criteria between $\Phi_{POD}$ and $\Phi$')
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

figure(6),clf
MAC(som_smd, real(Phi_x));
%     MAC(som_smd, real(Phi_x(:,1:size(som_smd,2))));
sgtitle(['Modal Assurance Criteria between $\Phi_{SCD}$ and $\Phi$; rr = ',num2str(rr)])
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])
%% Transformed Smooth Mode Decomposition - spatial dimension reduction
r = size(Y,2);
r = 100;
% Direct POD
[U, S, V] = svd([Y; delY]);
% Projection matrix
Vt = V(:,1:r);
% Data projection
Yt = Y*Vt;
delYt = delY*Vt;
% TSMD
[somt, sovt, spmt, soct, S1t, S2t, ~, ~] = sod(Yt', delYt');
% Project back
som = Vt*somt;

% Apply the modal assurance criterion to obtain the sorting indices
[~, ~, SIndx_smd] = MAC(som, real(Phi_x));
[~, ~, SIndx_pod] = MAC(pom, real(Phi_x));

% ------------------- VISUALIZATION OF ESTIMATED PHI --------------------
% Plot estimation results
figure(2),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    pom_n = normalize(pom(:,SIndx_pod(i)),'norm');
    som_n = normalize(som(:,SIndx_smd(i)),'norm');
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
    ylim([-0.05 0.05])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])

sgtitle('Mode Shape Estiamtion between POD and SCD')


% Plot the comparison between POD and SMD
nrows = 8;
figure(3),clf
for i = 1:nrows
subplot(nrows,2,2*(i-1)+1)
imagesc(vec2snap(som(:,SIndx_smd(i)),mm,nn)), mycolorbar()
pbaspect([b0 a0 1])
subplot(nrows,2,2*(i-1)+2)
plot(soct(:,SIndx_smd(i)))
pbaspect([b0 a0 1])
end

nrows = 8;
figure(4),clf
for i = 1:nrows
subplot(nrows,2,2*(i-1)+1)
imagesc(vec2snap(pom(:,SIndx_pod(i)),mm,nn))
pbaspect([b0 a0 1])
subplot(nrows,2,2*(i-1)+2)
plot(poc(:,SIndx_pod(i)))
pbaspect([b0 a0 1])
end

% ------------------------------------------------------------------------
%                 Temporally underdetermined problem 
% ------------------------------------------------------------------------
%% Temporal Dimensionality Reduction (TDR-3) --- TSVD-DEL based
% ------------------------------------------------------------------------
%        Transformed/Truncated Smooth Mode Decomposition (Temporal)
%      --------------------------------------------------------------
% 1. Temporal dimensionality reduction by obtaining an orthogonal
%    projector, POC of [X, delX], and truncated it to some degree
% 2. Data projection onto the common POM, Ut
% 3. SMD using the projected data sets [Xt, delXt]
% 4. Unfold the modes to its original dimensionality 
% ------------------------------------------------------------------------
PrintFlag = 0;
SpectralPlot = 1;
Min_rr_SMD = m/n;
% Set up the resampling rate 
rr = 1; % TO GENERATE THE ONE USED IN THE REPORT rr = 460;
var_Y_mean = mean(var(Y_raw));
SNR = 1000;
NoiseLevel = var_Y_mean/SNR;
% if rr < Min_rr_SMD
%     error(['The minimum size of resampling rate has to be larger than:',num2str(Min_rr_SMD)])
% end
rng(1)
% Data Resampling (down-sampling)
Y_raw_noise = Y_raw + NoiseLevel*randn(size(Y_raw));
Y = Y_raw_noise(1:rr:end, :); 
% Put data into a tensor
TY = Y2TY(Y, mm, nn);
[delY, ~] = TY2delY(TY, dx, dy);

% 1. Direct POD to concatenated matrices [Y, DY]
[poc, pov, pom] = svd(Y_raw, 'econ');

% DY = fs*diff(Y);
[Ut, St, Vt] = svd([Y,delY], 'econ');


% 2. TSVD based temporal dimensionality reduction
r = 10;
Q = Ut(:, 1:r);
Yr = Q'*Y;
delYr = Q'*delY; % Center finite difference

% 3. TSMD
[som_tsmd, sov_tsmd, spm_tsmd, sc_tsmd, S1_tsmd, S2_tsmd, U_tsmd, V_tsmd] = sod(Yr', delYr', ' ', 'False');

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
    plot(dir_pod*pom_n)
    hold on

    % SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS
    plot(dir_sod*som_n)
    plot(Phi_am_n, 'k--')
    axis tight

    ylabel(['Mode ',num2str(i)])

    if count == 1
        legend(['$\widehat\phi_\mathrm{pod}$'],['$\widehat\phi_\mathrm{smd}$'],['$\phi$'],'NumColumns',3, 'location','northoutside')
    end

    if count == 9 || count == 10
        xlabel('Noal Point')
    end
    count = count + 1;
    ylim([-0.05 0.05])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])

sgt = sgtitle(['Mode Shape Estiamtion between POD and TSMD-N; $rr = ',num2str(rr),';r = ',num2str(r),'$']);
sgt.FontSize = 12;

figure(3),clf,MAC(pom, real(Phi_x));
sgtitle('Modal Assurance Criteria between $\Phi_{POD}$ and $\Phi$')
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

subplot(121)
xlim([0.5 r + 0.5])
ylim([0.5 r + 0.5])
subplot(122)
xlim([0.5 r + 0.5])
ylim([0.5 r + 0.5])

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

% Plot estimation results (cooridnates)
figure(5),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    poc_n = normalize(poc(1:end,SIndx_pod(i)),'norm');
    sc_n = normalize(sc(:,SIndx_tsmd(i)),'norm');
    Eta_n = normalize(Eta(1:end,i),'norm');
    % Check the orientation of the modes
    dir_poc = sign(poc_n'*Eta_n);
    % Check the orientation of the modes
    dir_sc = sign(sc_n'*Eta_n);
    
    windsize = round(m/3);
    noverlap = round(windsize/2);
    [pxx_poc, fxx] = pwelch(dir_poc*poc_n, hanning(windsize), noverlap, [], fs);
    [pxx_sc, ~] = pwelch(dir_sc*sc_n, hanning(windsize), noverlap, [], fs);
    [pxx_eta, ~] = pwelch(Eta_n, hanning(windsize), noverlap, [], fs);
    if ~SpectralPlot
        plot([dir_poc*poc_n])
        hold on
        plot(dir_sc*sc_n)
        plot(Eta_n,'--k')
        if count == 9 || count == 10
            xlabel('Sample Time')
        end
        axis tight
    else
        plot(fxx, 10*log10(pxx_poc))
        hold on
        plot(fxx, 10*log10(pxx_sc))
        plot(fxx, 10*log10(pxx_eta), 'k--')
        if count == 9 || count == 10
            xlabel('Frequency (Hz)')
        end
        xlim([0 5])
    end
    % SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS

    ylabel(['Mode ',num2str(i)])

    if count == 1
        legend(['$\widehat\eta_\mathrm{pod}$'],['$\widehat\eta_\mathrm{tsmd}$'],['$\eta$'],'NumColumns',3, 'location','northoutside')
    end

    count = count + 1;
    %     xlim([0 5])
    %     ylim([-0.02 0.02])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])
sgt = sgtitle(['Modal Coordinate Estiamtion between POD and TSMD-N; $rr = ',num2str(rr),';r = ',num2str(r),'$']);
sgt.FontSize = 12;


% Plot the comparison between POD and SMD
nrows = 8;
figure(3),clf
for i = 1:nrows
subplot(nrows,2,2*(i-1)+1)
imagesc(vec2snap(som_tsmd(:,SIndx_tsmd(i)),mm,nn)), mycolorbar()
pbaspect([b0 a0 1])
subplot(nrows,2,2*(i-1)+2)
plot(sc(:,SIndx_tsmd(i)))
pbaspect([b0 a0 1])
end

nrows = 8;
figure(4),clf
for i = 1:nrows
subplot(nrows,2,2*(i-1)+1)
imagesc(vec2snap(pom(:,SIndx_pod(i)),mm,nn))
pbaspect([b0 a0 1])
subplot(nrows,2,2*(i-1)+2)
plot(poc(:,SIndx_pod(i)))
pbaspect([b0 a0 1])
end

if PrintFlag == 1
    figure(2), print('POD_TSMDN_phi_comparison_pp_beam.png','-dpng','-r600')
    figure(4), print(['MAC_TSMDN_rr_',num2str(rr),'_r_',num2str(r),'.png'],'-dpng','-r600')
    figure(5), print('POD_TSMDN_eta_comparison_pp_beam.png','-dpng','-r600')
end
%% Truncated Smooth Coordinate Decomposition
% ------------------------------------------------------------------------
%              TRUNCATED SMOOTH COORDINATE DECOMPOSITION
%            ---------------------------------------------
% 1. Spatial demensionality reduction by obtaining an orthgonal projector,
%    the common POM using the TSVD of [X; delX]
% 2. Data projection onto the common POM, Vt
% 3. SCD using the projected data sets [Xt, DXt]
% 4. Unfold the modes to its original dimensionality
% ------------------------------------------------------------------------
PrintFlag = 0;
SpectralPlot = 1;
Difference = 1;
Min_rr_SMD = m/n;
% Set up the resampling rate 
rr = 1; % TO GENERATE THE ONE USED IN THE REPORT rr = 460;
% var_Y_mean = mean(var(Y_raw));
SNR = 100;
% NoiseLevel = var_Y_mean/SNR;

% Data Resampling (down-sampling)
rng(1)
Y_raw_noise = awgn(Y_raw, SNR, 'measured',[],'dB');
Y = Y_raw_noise(1:rr:end, :); 

% Windowing Parameter
windsize = size(Y,1)/2;
noverlap = windsize/4;
wind = hanning(windsize);

% Put data into a tensor
TY = Y2TY(Y, mm, nn);
[delY, ~] = TY2delY(TY, dx, dy);


% ------------------- VISUALIZATION OF RESAMPLED DATA --------------------
figure(1),clf
imagesc(Y)
view(270,-90)
pbaspect([1 2.5 1])
ylabel('Time Samples')
xlabel('Spatial Samples')
set(gcf,'papersize', [6 2.5])
set(gcf,'paperposition', [0 0 6 2.5])
mycolormap = mycolorbar('Viridis');
colormap(mycolormap)

DY = fs*diff(Y);
% Spatial Dimensionality Reduction through TSVD
r = size(Y,1) - 10;
r = 10;
% Direct POD
[U, S, V] = svd([Y; delY]);
% Projection matrix
Vt = V(:,1:r);
% Data projection
Yt = Y*Vt;
DYt = DY*Vt;
[soc_tscd, sov_tscd, spm_tscd, sm_tscd, S1, S2, U1, U2] = sod(Yt, DYt, ' ', 'false');

sm = Vt*sm_tscd;

[poc, pov, pom] = svd(Y, 'econ');



% Apply the modal assurance criterion to obtain the sorting indices
[~, MIndx_scd, SIndx_scd] = MAC(sm, real(Phi_x));
[~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));

% ------------------- VISUALIZATION OF ESTIMATED PHI --------------------
% Plot estimation results
figure(2),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    pom_n = normalize(pom(:,SIndx_pod(i)),'norm');
    som_n = normalize(sm(:,SIndx_scd(i)),'norm');
    Phi_am_n = normalize(Phi_am(:,i),'norm');
    % Check the orientation of the modes
    dir_pod = sign(pom_n'*Phi_am_n);
    % Check the orientation of the modes
    dir_sod = sign(som_n'*Phi_am_n);

    pom_n_dir = dir_pod*pom_n;
    som_n_dir = dir_sod*som_n;
    d_pom_n_dir = Phi_am_n - pom_n_dir;
    d_som_n_dir = Phi_am_n - som_n_dir;
    
    if ~Difference
        plot(pom_n_dir)
        hold on
        plot(som_n_dir)

        plot(Phi_am_n, 'k--')
        axis tight
    else
        plot(d_pom_n_dir)
        hold on
        plot(d_som_n_dir)

        plot(zeros(length(d_pom_n_dir),1), 'k--')
        axis tight
    end
    ylabel(['Mode ',num2str(i)])
    
    if count == 1
        legend(['$\widehat\phi_\mathrm{POD}$'],['$\widehat\phi_\mathrm{TSCD}$'],['$\widehat\phi$'],'NumColumns',3, 'location','northoutside')
    end

    if count == 9 || count == 10
        xlabel('Noal Point')
    end
    count = count + 1;
%     ylim([-0.05 0.05])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])
if ~Difference
    sgtitle('Mode Shape Estiamtion between POD and TSCD')
else
    sgtitle('Mode Shape Estiamtion Error between POD and TSCD')
end

% --------------- VISUALIZATION OF THE ESTIMATED PHI --------------------

figure(3),clf
count = 1;
nrows = 10;
ncols = 3;
for i = 1:nrows
    % Get the normalized modes (according to some normalization scheme)
    pom_n = normalize(pom(:,SIndx_pod(i)),'norm');
    som_n = normalize(sm(:,SIndx_scd(i)),'norm');
    Phi_am_n = normalize(Phi_am(:,i),'norm');
    % Check the orientation of the modes
    dir_pod = sign(pom_n'*Phi_am_n);
    % Check the orientation of the modes
    dir_sod = sign(som_n'*Phi_am_n);
    % SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS

    % First column of plot
    subplot(nrows,ncols,ncols*(i-1) + 1)
    imagesc(reshape(Phi_am_n,mm,nn))
    xticks([]);
    yticks([]);
    ylabel(['Mode ',num2str(i)])
    if i == 1
        title('Truth - $\Phi_i$')
    end

    mycolorbar()
    % Second column of plot
    subplot(nrows,ncols,ncols*(i-1) + 2)
    imagesc(reshape(dir_pod*pom_n,mm,nn))
    xticks([]);
    yticks([]);
    if i == 1
        title('POM - $\Phi^\mathrm{POD}_{i}$')
    end

    % Third column of plot
    subplot(nrows,ncols,ncols*(i-1) + 3)
    imagesc(reshape(dir_sod*som_n,mm,nn))
    xticks([]);
    yticks([]);
    if i == 1
        title('SM - $\Phi^\mathrm{TSCD}_{i}$')
    end
    axis tight

    count = count + 1;
    %     ylim([-0.05 0.05])
end
set(gcf, 'papersize', [10 14])
set(gcf, 'paperposition', [0 0 10 14])

sgtitle('Mode Shape Estiamtion between POD and TSCD')

% ------------------ MAC FOR POD MODES -----------------------------------
figure(4),clf,MAC(pom, real(Phi_x));
sgtitle('Modal Assurance Criteria between $\Phi_\mathrm{POD}$ and $\Phi$')
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

subplot(121)
xlim([0.5 10 + 0.5])
ylim([0.5 10 + 0.5])
subplot(122)
xlim([0.5 10 + 0.5])
ylim([0.5 10 + 0.5])

% ------------------ MAC FOR TSCD MODES ----------------------------------
figure(5),clf
MAC(sm, real(Phi_x));
sgtitle(['Modal Assurance Criteria between $\Phi_\mathrm{TSCD}$ and $\Phi$; r = ',num2str(r)])
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

subplot(121)
xlim([0.5 10 + 0.5])
ylim([0.5 10 + 0.5])
subplot(122)
xlim([0.5 10 + 0.5])
ylim([0.5 10 + 0.5])

% ------------------- VISUALIZATION OF ESTIMATED ETA --------------------
figure(6),clf
count = 1;
strt = 1;
for i = strt:strt+9
    subplot(5,2,count)
    % Get the normalized modes (according to some normalization scheme)
    poc_n = normalize(poc(1:end,SIndx_pod(i)),'norm');
    sc_n = normalize(soc_tscd(:,SIndx_scd(i)),'norm');
    Eta_n = normalize(Eta(1:end,i),'norm');
    % Check the orientation of the modes
    dir_poc = sign(poc_n'*Eta_n);
    % Check the orientation of the modes
    dir_sc = sign(sc_n'*Eta_n);
    
    windsize = round(m/4);
    noverlap = round(windsize/2);
    [pxx_poc, fxx] = pwelch(dir_poc*poc_n, wind, noverlap, [], fs);
    [pxx_sc, ~] = pwelch(dir_sc*sc_n, wind, noverlap, [], fs);
    [pxx_eta, ~] = pwelch(Eta_n, wind, noverlap, [], fs);
    if ~SpectralPlot
        plot([dir_poc*poc_n])
        hold on
        plot(dir_sc*sc_n)
        plot(Eta_n,'--k')
        if count == 9 || count == 10
            xlabel('Sample Time')
        end
        axis tight
    else
        plot(fxx, 10*log10(pxx_poc))
        hold on
        plot(fxx, 10*log10(pxx_sc))
        plot(fxx, 10*log10(pxx_eta), 'k--')
        if count == 9 || count == 10
            xlabel('Frequency (Hz)')
        end
        xlim([0 5])
    end
    % SCD - IF SCD IS DIRECTLY USED, IT IS GOING TO FAIL DUE TO ILL-POSENESS

    ylabel(['Mode ',num2str(i)])

    if count == 1
        legend(['$\widehat{q}_\mathrm{POD}$'],['$\widehat{q}_\mathrm{TSCD}$'],['$q$'],'NumColumns',3, 'location','northoutside')
    end

    count = count + 1;
    %     xlim([0 5])
    %     ylim([-0.02 0.02])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])
sgt = sgtitle(['Modal Coordinate Estiamtion between POD and TSCD; $r = ',num2str(r),'$']);
sgt.FontSize = 12;

% ------------------- PRINT RESULTS --------------------
if PrintFlag == 1
    figure(1), print(['Y_resampled_rr_',num2str(rr),'_',CurrentCase,'.png'],'-dpng','-r600')
    figure(2), print(['POD_TSCD_phi_comparison_ffm_',CurrentCase,'.png'],'-dpng','-r600')
    figure(3), print(['POD_TSCD_Phi_comparison_ffm_',CurrentCase,'.png'],'-dpng','-r600')
    figure(4), print(['MAC_POD_r_',num2str(r),'_',CurrentCase,'.png'],'-dpng','-r600')
    figure(5), print(['MAC_TSCD_r_',num2str(r),'_',CurrentCase,'.png'],'-dpng','-r600')
    figure(6), print(['POD_TSCD_eta_comparison_pp_beam_',CurrentCase,'.png'],'-dpng','-r600')
end