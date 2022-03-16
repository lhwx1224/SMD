
%% UNDERDETERMINED SCD
clear U S V
Y = y_modal;
DY = fs*diff(Y);                          % Temporal Differentiation
[U, S, V] = svd([Y; DY], 'econ'); 
[Up, Sp, Vp] = svd(Y, 'econ'); 

% Decide where to truncate the rank
figure(2),clf
plot(log10(diag(S)))
% Spatial Truncation
r = 180;
Vr = V(:,1:r);  % New basis for the data
Vrp = Vp(:,1:r);
xs = round(linspace(1, length(Phi_am), r));
Yt = Y*Vr;      % Project the data onto the new basis
DYt = DY*Vr;    % Project the derivative onto the same basis

% SVD to the raw data
[poct, povt, pomt] = svd(Yt, 'econ');
% SCD
addpath('L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox')
[soct, sovt, spmt, somt, S1t, S2t, ~, ~] = sod(Yt, DYt);

som = Vr*somt;
% pom = Vr*pomt;
pom = Vrp*pomt;
pom = fliplr(pom);
% 
% subplot(211)
% plot(soct(:,3))
% subplot(212)
% plot(som(:,3))

figure(3),clf
count = 1;
for i = 1:10
subplot(5,2,count)
% Get the normalized modes (according to some normalization scheme)
pom_n = normalize(pom(:,i),'norm');
som_n = normalize(som(:,i),'norm');
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
legend(['$\widehat\phi_\mathrm{pod}$'],['$\widehat\phi_\mathrm{sod}$'],['$\widehat\phi$'],'NumColumns',3, 'location','northoutside')
end

if count == 9 || i == 10
    xlabel('Noal Point')
end
count = count + 1;
end
sgtitle('Mode Shape Estiamtion between POD and SOD (RSCD)')


set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])

MAC(som, real(Phi_x))
sgtitle('Modal Assurance Criteria between $\Phi_{RSCD}$ and $\Phi$')
MAC(pom, real(Phi_x))
sgtitle('Modal Assurance Criteria between $\Phi_{POD}$ and $\Phi$')
%% SMOOTH MODE DECOMPOSITION
% ------------------------------------------------------------------------
% SMD of a spatial-temporally underdetermined problem cannot resolve unique
% solution of the system
% ------------------------------------------------------------------------
% Set up the environment
addpath('L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox')
% Set up the resampling rate 
rr = 460;
% Data Resampling (down-sampling)
Y = y_modal(1:rr:end, :); 
% Y = Y + 0.0001*randn(size(Y));    % Add noise (optional)
% Apply finite difference method to the data
delY = GenFiniteDiff(Y', dx, 'c2')'; % Center difference with padding

% Direct POD
[poc, pov, pom] = svd(Y, 'econ');

% DIRECT SMD
tic
[som_smd, sov_smd, spm_smd, soc_smd, S1_smd, S2_smd, U_smd, V_smd] = sod(Y', delY');
toc

% Apply the modal assurance criterion to obtain the sorting indices
[~, MIndx_smd, SIndx_smd] = MAC(som_smd, real(Phi_x));
[~, MIndx_pod, SIndx_pod] = MAC(pom, real(Phi_x));

% Plot estimation results
figure(4),clf
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
ylim([-0.07 0.07])
end
set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])
print('POD_SMD_comparison_pp_beam.png','-dpng','-r600')

sgtitle('Mode Shape Estiamtion between POD and SOD (SMD)')
figure(5),clf,MAC(pom, real(Phi_x));
sgtitle('Modal Assurance Criteria between $\Phi_{POD}$ and $\Phi$')
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

figure(6),clf,MAC(som_smd, real(Phi_x(:,1:size(som_smd,2))));
sgtitle('Modal Assurance Criteria between $\Phi_{SMD}$ and $\Phi$')
set(gcf, 'papersize', [6 3])
set(gcf, 'paperposition', [0 0 6 3])

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
    Y = y_modal(1:50+j:end, :); % 188， 210， 267, 310, 400, 464, 479
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
ylabel('Cumulative Error - $\sum_{i = 1}^{50}(\hat\phi_i - \phi_i)^2$')
grid on
set(gcf,'papersize',[6 2.5])
set(gcf,'paperposition',[0 0 6 2.5])
%% TSMD
% SMD
[poc, pov, pom] = svd(Y,'econ');
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

set(gcf, 'papersize', [6 6])
set(gcf, 'paperposition', [0 0 6 6])
print('POD_RSMD_comparison_pp_beam.png','-dpng','-r600')

% sgtitle('Mode Shape Estiamtion between POD and SOD (RSMD)')
% Shrink the size in the spatial domain
figure(6),clf
MAC(som_smd, real(Phi_x))
sgtitle('Modal Assurance Criteria between $\Phi_{RSMD}$ and $\Phi$')
figure(7),clf
MAC(pom, real(Phi_x))
sgtitle('Modal Assurance Criteria between $\Phi_{POD}$ and $\Phi$')