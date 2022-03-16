%% SOD MAIN DATA POST-PROCESSING FILE
set(0,'DefaultFigureWindowStyle','docked')
DefaultRootPath = 'L:\My Drive\Graduate study\Research\Projects\';
RootPath = DefaultRootPath;
%% Load the video into the matrix
close all
addpath([RootPath,'Output-only Modal Analysis Toolbox'])
Y_raw = y_modal;
delY_raw = dely_modal;
%% ON-DEMAND SAVE DATA 
save(['Y_raw.mat'],'Y_raw','-v7.3')

%% Preview of the video             
StartFrame = 1;
% EndFrame = 200;
dFrame = 1;
for i = StartFrame:dFrame:dFrame*fs*2
   u1 = reshape(Y(i,:), Nx, Ny);
   imagesc(u1),colormap('bone')
   pbaspect([Ny/Nx 1 1])
   title(['Frame:', num2str(i), ' out of ',num2str(dFrame*fs*2), '; Display Rate:', num2str(dFrame),'.'])
   pause(0.0001)
end
%% Proper Orthogonal Decomposition (Full SVD)
Y = Y_raw;
windsize = round(size(Y,1)/2);
noverlap = windsize/2;

tic
[U, S, V] = svd(Y - mean(Y), 'econ');
toc

i = 1;
figure(1),clf
for k = 1:3
subplot(3,3,3*(k-1)+1)
imagesc(reshape(V(:,i), Nx, Ny))
pbaspect([b0/a0 1 1])
mycolorbar();
title(['POD-Mode#',num2str(i)])
% colorbar()
subplot(3,3,3*(k-1)+2)
plot(U(:,i),'-')
pbaspect([b0/a0 1 1])
title(['POD-Coordinate#',num2str(i)])
xlabel('Sample')
ylabel('Magnitude')
axis tight
subplot(3,3,3*(k-1)+3)
[pxx_pod, freq_pod] = pwelch(U(:,i),boxcar(round(windsize*size(Y,1)/size(Y_raw,1))),round(noverlap*size(Y,1)/size(Y_raw,1)),[],fs);
plot(freq_pod, 10*log10(pxx_pod))
[pxx_pod_max, indx] = max(pxx_pod);
hold on
plot(freq_pod(indx), 10*log10(pxx_pod_max),'o')
pbaspect([b0/a0 1 1])
grid on
title(['PSD of POD-Coordinate#',num2str(i)])
xlabel('Frequency (Hz)')
ylabel('PSD')
axis tight
i = i + 1;
end
figure(1)
set(gcf,'PaperSize',[20 12])
set(gcf,'PaperPosition',[0 0 20 12])

%% Truncated SVD
r = 5;                              % Order of truncation - r
Yr = Y*V(:,1:r);                    % Project the data onto the first r POD modes
[Ur, Sr, Vr] = svd(Yr, 'econ');     % TSVD
i = 2;                             
figure(1),clf
subplot(311)
imshow(reshape(V(:,1:r)*Vr(:,i), Nx, Ny), [])
pbaspect([Ny/Nx 1 1])
title(['POD-Mode#',num2str(i)])
% colorbar()
subplot(312)
plot(Ur(1:fs,i),'-')
pbaspect([Ny/Nx 1 1])
title(['POD-Coordinate#',num2str(i)])
xlabel('Sample')
ylabel('Magnitude')

subplot(313)
[pxx_pod, freq_pod] = pwelch(Ur(:,i),boxcar(windsize),noverlap,[],fs);
plot(freq_pod, 10*log10(pxx_pod))
[pxx_pod_max, indx] = max(pxx_pod);
hold on
plot(freq_pod(indx), 10*log10(pxx_pod_max),'o')
pbaspect([Ny/Nx 1 1])
grid on
title(['PSD of POD-Coordinate#',num2str(i)])
xlabel('Frequency (Hz)')
ylabel('PSD')
%% SINGULAR VALUE DECOMPOSITION 
% u1 = reshape(V(:,2),Nx, Ny);
% imagesc(u1),colormap('gray')
% pbaspect([Ny/Nx 1 1])
%% POD COMPONENT-WISE VIEW
Nt = 200;
Yrk = zeros(size(Y)); 
r = 1;
for k = r
    Yrk = Yrk + U(:,k)*S(k,k)*V(:,k)';
end
for i = 1:1:Nt
   subplot(211)
   Yr1_temp = Yrk(i,:);
   Yr1_temp = reshape(Yr1_temp, Nx, Ny);
   imshow(Yr1_temp,[0 255]),colormap('bone')
   pbaspect([Ny/Nx 1 1])
   title(['POD Mode #',num2str(r)])
   pause(0.001)
   sp2 = subplot(212);
   cla(sp2)
   plot(1:length(U(:,r)), U(:,r))
   hold on
   plot(i, U(i,r),'ro')
   pbaspect([Ny/Nx 1 1])
   xlim([0 Nt])
end
%% POD LOW-RANK Approximation
Yrk = zeros(size(Y));
r = 10;
for k = 1:r
    Yrk = Yrk + U(:,k)*S(k,k)*V(:,k)';
end
clf
for i = 1:1:Nt
   Yr1_temp = Yrk(i,:);
   Yr1_temp = reshape(Yr1_temp, Nx, Ny);
   imagesc(Yr1_temp),colormap('bone')
   pbaspect([Ny/Nx 1 1])
   title(['POD Rank-',num2str(r-1),' Approximation'])
   pause(0.001)
end

windsize = 2^nextpow2(size(Y,1)/2);

noverlap = windsize/2;

figure(2),clf
pxxm = zeros(10,1);
fm = zeros(10,1); 
for i = 1:50
[pxx, f] = pwelch(U(:,i),boxcar(windsize),noverlap,[],fs);
[pxxm(i), pxxmi(i)] = max(pxx);
fm(i) = f(pxxmi(i));
plot(fm(i),10*log10(pxxm(i)),'k.')
hold on
end
[pxx, f] = pwelch(U(:,2),boxcar(windsize),noverlap,[],100);
plot(f, 10*log10(pxx))
%% EXTRACT THE MODE SHAPE FROM THE IDENTIFICATION DATA

set(gcf, 'renderer','opengl')
figure(1),clf
u1 = reshape(V(:,2), Nx, Ny);
for i = 1:1:size(u1,2)
    plt = u1(1:1:end,i);
    [pks, pindx]= max(plt);
    plot3(i*ones(length(plt),1),linspace(0,length(plt),length(plt)),plt,'k.','MarkerSize',0.5)
    hold on
%     plot3(i*ones(length(diff(plt)),1),linspace(0,length(diff(plt)),length(diff(plt))),diff(plt),'b*','MarkerSize',6)
    plot3(i*ones(length(pks),1),pindx,pks,'ro')
end
pbaspect([Ny/Nx 1 1])
view([0 90])
%% Plot only extracted mode shape
figure(1),clf
pks = zeros(size(u1,2),1);
pindx = zeros(size(u1,2),1);
for i = 1:1:size(u1,2)
    plt = u1(1:2:end, i);
    [pks(i), pindx(i)]= max(plt);
    plot3(i*ones(length(pks(i)),1),pindx(i),pks(i),'ro')
    hold on    
end
pbaspect([Ny/Nx 1 1])
view([0 90])
% xlim([41 1360])

%% Curve fitting to the mode shape
height = 1:length(pindx(41:1360));
width = pindx(41:1360);

nheight = FirstLast_normalize(height);
nwidth = flipud(FirstLast_normalize(width));
nwidth = nwidth - nwidth(1);

plot(nheight, nwidth ,'o')
pbaspect([Ny/Nx 1 1])

[xData, yData] = prepareCurveData( nheight, nwidth );

% Set up fittype and options.
ft = fittype( '(sin(b*x) - sinh(b*x) - ((sin(b) + sinh(b))/(cos(b) + cosh(b)))*(cos(b*x) - cosh(b*x)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.655740699156587];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'nwidth vs. nheight', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'nheight', 'Interpreter', 'none' );
ylabel( 'nwidth', 'Interpreter', 'none' );
grid on

%% UNDERDETERMINED SCD
% ------------------------------------------------------------------------
% For underdetermined problem, use smooth coordinate decomposition directly
% with truncated SVD
% ------------------------------------------------------------------------
clear U S V
DY = fs*diff(Y);                          % Temporal Differentiation
imagesc(DY)
[U, S, V] = svd([Y; DY], 'econ'); 

figure(1),clf
subplot(211)
plot(diag(S)), set(gca,'yscale','log')
pbaspect([2 1 1])
xlabel('Index of Singular Values')
ylabel('Singular Values')
grid on
subplot(212)
plot(diag(S)), set(gca,'yscale','log')
pbaspect([2 1 1])
xlim([0 100])
xlabel('Index of Singular Values')
ylabel('Singular Values')
grid on
figure(1), set(gcf,'papersize',[6 6]),set(gcf,'paperposition',[0 0 6 6])

% Spatial Truncation
r = 100;
Vr = V(:,1:r);  % New basis for the data

Yt = Y*Vr;      % Project the data onto the new basis
DYt = DY*Vr;    % Project the derivative onto the same basis
% SCD
addpath('L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox')
[soct, sovt, spmt, somt, S1t, S2t, Ut, Vt] = sod(Yt, DYt);

som = Vr*somt;

i = 1;
subplot(311)
imagesc(reshape(som(:,i), Nx, Ny)), colormap('gray')
pbaspect([Ny/Nx 1 1])
xlabel('x')
ylabel('y(x)')
title(['SCD-Mode#',num2str(i)])
subplot(312)
plot(soct(:,i))
title(['SDC-Coordinate#',num2str(i)])
xlabel('Samples')
ylabel('Magnitude')
pbaspect([Ny/Nx 1 1])
axis tight
subplot(313)
[pxx, freq] = pwelch(soct(:,i), boxcar(windsize), noverlap, [], fs);
plot(freq, 10*log10(pxx))
title(['Pxx of SDC-Coordinate#',num2str(i)])
xlabel('Frequency (Hz)')
ylabel('PSD')
pbaspect([Ny/Nx 1 1])
figure(1),set(gcf,'PaperSize',[6 6*1.5]),set(gcf,'PaperPosition',[0 0 6 6*1.5])

%% UNDERDETERMINED SCD: Low-rank Approximation Animation
%
%
%
clf
Ysk = zeros(size(Y));
r = 10;
for k = 1:r
    Ysk = Ysk + soct(:,k)*som(:,k)';
end
figure(1)
colorbar()
for i = 1:1:100
   Ys1_temp = Ysk(i,:);
   Ys1_temp = reshape(Ys1_temp, Nx, Ny);
   imshow(Ys1_temp,[-255 255])
   title(['SCD Rank-',num2str(r),' Approximation'])
   pbaspect([Ny/Nx 1 1])
   pause(0.001)
end

%% SMOOTH MODE DECOMPOSITION TO THE VIDEO DATA SETS
% ----------------------------------------------------------------
% Full SMD
% ----------------------------------------------------------------
rr = 60;
EndFrame = size(Y_raw,1);
Y = Y_raw(1:rr:EndFrame,:);
delY = delY_raw(1:rr:EndFrame,:);
% Y = Y_raw(unique(randi([1 size(Y_raw,1)], 600, 1)), :);
addpath('L:\My Drive\Graduate study\Research\My_paper\Journal\MSSP\LNN_NNM_ID\Code\')
% Preprocess the video data
% delY = zeros(size(Y));    % Prepare for the differentiated data
% Finite Difference in the horizontal direction
% for i = 1:size(Y, 1)
%     Y_temp = reshape(Y(i, :), Nx, Ny);            % assign temporary reshaped matrix
%     delY_temp = GenFiniteDiff(Y_temp', 1, 'c2')'; % central difference 
%     delY(i,:) = reshape(delY_temp, (Nx)*(Ny), 1)';% assign the differentiated-and-reshaped data
% end
% Visualization to the differentiated field
% for i = 1:1:20
%     % Plot the original video snapshot
%     subplot(211)
%     u1 = reshape(Y(i,:), Nx, Ny);
%     imagesc(u1),colormap('bone')
%     pbaspect([Ny/Nx 1 1])
%     title('The Original Snapshot')
%     % Plot the differentiated video snapshot 
%     subplot(212)
%     u1 = reshape(delY(i,:), Nx, Ny);
%     imagesc(u1),colormap('bone')
%     pbaspect([Ny/Nx 1 1])
%     title('Spatially Differentiated Snapshot')
%     pause(0.01)
% end
tic
[somt, sovt, spmt, soct, S1t, S2t, Ut, Vt] = sod(Y', delY');
toc

figure(1),clf
subplot(211)
plot(sovt), set(gca,'yscale','log')
pbaspect([2 1 1])
xlabel('Index of Singular Values')
ylabel('SOV')
title('TSMD - Smooth Orthogonal Values')
grid on
subplot(212)
plot(sovt), set(gca,'yscale','log')
pbaspect([2 1 1])
xlim([0 50])
xlabel('Index of SOVs')
ylabel('Singular Values')
title('TSMD - Smooth Orthogonal Values (close-up)')
grid on
figure(1), set(gcf,'papersize',[6 6]),set(gcf,'paperposition',[0 0 6 6])

% i = 1;
% figure(2),clf
% subplot(311)
% som_SMD = reshape(somt(:,i), Nx, Ny);
% imshow(som_SMD,[]), colormap('gray')
% pbaspect([Ny/Nx 1 1])
% xlabel('x')
% ylabel('y(x)')
% title(['SMD-Mode#',num2str(i)])
% subplot(312)
% plot(soct(:,i))
% title(['SMD-Coordinate#',num2str(i)])
% xlabel('Samples')
% ylabel('Magnitude')
% pbaspect([Ny/Nx 1 1])
% subplot(313)
% [pxx, freq] = pwelch(soct(:,i), hanning(windsize), noverlap, [], 100);
% plot(freq, 10*log10(pxx))
% title(['Pxx of SMD-Coordinate#',num2str(i)])
% xlabel('Frequency (Hz)')
% ylabel('PSD')
% pbaspect([Ny/Nx 1 1])
% figure(2),set(gcf,'PaperSize',[6 6*1.5]),set(gcf,'PaperPosition',[0 0 6 6*1.5])

figure(2),clf
i = 1;
for k = 1:3
subplot(3,3,3*(k-1)+1)
som_SMD = reshape(somt(:,i), Nx, Ny);
imshow(som_SMD,[]), colormap('gray')
pbaspect([Ny/Nx 1 1])
% xlabel('x')
% ylabel('y(x)')
title(['SMD-Mode#',num2str(i)])
subplot(3,3,3*(k-1)+2)
plot(soct(:,i))
title(['SDC-Coordinate#',num2str(i)])
xlabel('Samples')
ylabel('Magnitude')
pbaspect([Ny/Nx 1 1])
axis tight
subplot(3,3,3*(k-1)+3)
[pxx, freq] = pwelch(soct(:,i), boxcar(round(size(Y,1)/2)), round(round(size(Y,1)/2)/2), [], fs/rr);
plot(freq, 10*log10(pxx))
title(['Pxx of SDC-Coordinate#',num2str(i)])
xlabel('Frequency (Hz)')
ylabel('PSD')
pbaspect([Ny/Nx 1 1])
i = i + 1;
end
figure(2),set(gcf,'PaperSize',[20 12]),set(gcf,'PaperPosition',[0 0 20 12])


% plot(soct(:,1))

%% ----------------------------------------------------------------
% TSMD - Reduce Spatial Dimension to Save Time
% ----------------------------------------------------------------
figure(1),clf
subplot(211)
plot(sovt), set(gca,'yscale','log')
pbaspect([2 1 1])
xlabel('Index of SOV')
ylabel('SOV')
title('SMD - Smooth Orthogonal Values')
grid on
subplot(212)
plot(sovt), set(gca,'yscale','log')
pbaspect([2 1 1])
xlim([0 50])
xlabel('Index of SOV')
ylabel('SOV')
title('SMD - Smooth Orthogonal Values (close-up)')
grid on
figure(1), set(gcf,'papersize',[6 6]),set(gcf,'paperposition',[0 0 6 6])

% RSVD
l = 110;                           % Sampling width (rank of the reduced input space)
O = randn(size([Y; delY], 1), l);   % Random sampling matrix O \in \mathbb{R}^{m_1 + m_2 \times l}
Yrp = [Y; delY]'*O;                 % Ramdomly sample the input space (rows of Y = [Y; delY])
[Q, R] = qr(Yrp,0);                 % Obtain the projection matrix (n-by-l)
Yt = Y*Q;                           % Project the data down to l-dimensional subspace
delYt = delY*Q;                     % Project the differentiated data down to l-dimensional subspace

% SMD
addpath('L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox')
tic
[somt, sovt, spmt, soct, S1t, S2t, Ut, Vt] = sod(Yt', delYt');
som_RSMD = Q*somt;
toc

% -------------- VISUALIZATION OF THE MODE -------------------------------
% Assign index of the resolved mode
i = 1;
subplot(311)
som_SMD = reshape(som_RSMD(:,i), Nx, Ny);
imshow(som_SMD,[]), colormap('gray')
pbaspect([Ny/Nx 1 1])
xlabel('x')
ylabel('y(x)')
title(['TSMD-Mode#',num2str(i)])
subplot(312)
plot(soct(:,i))
title(['TSMD-Coordinate#',num2str(i)])
xlabel('Samples')
ylabel('Magnitude')
pbaspect([Ny/Nx 1 1])
subplot(313)
[pxx, freq] = pwelch(soct(:,i), hanning(windsize), noverlap, [], 100);
plot(freq, 10*log10(pxx))
title(['Pxx of TSMD-Coordinate#',num2str(i)])
xlabel('Frequency (Hz)')
ylabel('PSD')
pbaspect([Ny/Nx 1 1])
figure(1),set(gcf,'PaperSize',[6 6*1.5]),set(gcf,'PaperPosition',[0 0 6 6*1.5])

%%
tic
[soc, sov, spm, som, Sigma1, Sigma2, U1, U2] = sod(Y, delY, ' ', 'False');
toc

plot(som(:,4))

for i = 4
    u1 = reshape(soc(:,i), Ny, Nx);
    imagesc(u1),colormap('bone')
    pbaspect([Nx/Ny 1 1])
    pause(0.01)
end

Ypk = zeros(size(Y));
for k = 2:4
    Ypk = Ypk + U(:,k)*S(k,k)*V(:,k)';
end
for i = 1:1:100
   Yp1_temp = Ypk(:,i);
   Yp1_temp = reshape(Yp1_temp, Ny, Nx);
   imagesc(Yp1_temp),colormap('bone')
   pbaspect([Nx/Ny 1 1])
   pause(0.001)
end

Ysk = zeros(size(Y));
for k = 1:3
    Ysk = Ysk + soc(:,k)*som(:,k)';
end
for i = 1:1:100
   Ys1_temp = Ysk(:,i);
   Ys1_temp = reshape(Ys1_temp, Ny, Nx);
   imagesc(Ys1_temp),colormap('bone')
   pbaspect([Nx/Ny 1 1])
   pause(0.001)
end

pwelch(V(:,2), boxcar(256), 128, [], 240)
hold on
pwelch(som(:,4) , boxcar(256), 128, [], 240)
%%
Ypk = zeros(size(Y));
Ysk = zeros(size(Y));

for k = 1:10
    Ypk = Ypk + U(:,k)*S(k,k)*V(:,k)';
end
for k = 1:10
    Ysk = Ysk + Sigma1(k)*U1(:,k)*som(:,k)';
end
Nframe = 100;
for i = 1:1:Nframe
    Yp1_temp = Ypk(:,i);
    Yp1_temp = reshape(Yp1_temp, Ny, Nx);
    Ys1_temp = Ysk(:,i);
    Ys1_temp = reshape(Ys1_temp, Ny, Nx);
    subplot(121)
    imagesc(Yp1_temp),colormap('gray')
    pbaspect([Nx/Ny 1 1])
    title('POD')
    subplot(122)
    imagesc(Ys1_temp),colormap('gray')
    title('SMD')
    pbaspect([Nx/Ny 1 1])
    sgtitle(['Frame: ',num2str(i),' of ', num2str(Nframe),'.'])
    pause(0.00001)
end
%%