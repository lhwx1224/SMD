function [yamp_test, yamp_FRF] = AccAmplifier(y, y_test, CHref, fs, fc, Method)
% AccAmplifier amplifies the modal test data 'y_test' based on a set of
% calibration data set 'y' (which may be calibrated first). A reference
% channel 'CHref' is used to designate the targeted channel signal based 
% off of which the amplifier is designed. Sampling parameters, sampling
% frequency 'fs' and cut-off frequency 'fc', are required for generating
% plots. 'Method' is used to choose from either Welch's method or FFT
% method for estimating the frequency response (FRF) functions.
% 
% Syntax: [yamp, yamp_FRF] = AccAmplifier(y, y_test, CHref, fs, fc, Method);
% 
% Input: y - 2d array, calibration data set
%        y_test - 2d array, modal test data set
%        CHref - int, reference channel index
%        fs - double, sampling frequency
%        fc - double, cut-off frequency
%        Method - String, (1) 'TF', transfer function estimate
%                         (2) 'FFT', fast Fourier transform estimate
%
% Output: yamp_test - 2d array, amplified modal test data
%         yamp_FRF - 2d array, frequency response function estimate for
%                    amplification 
%
% Required functions: DataInspection.m 
% Hewenxuan Li, Created on April 23, 2021
% Last Modified: 
% October 1, 2021: Added FFT-based amplitude amplification
% October 4, 2021: fully functional to directly incorporate both
% calibration data set and modal testing data set.

% Add toolbox path
addpath('L:\My Drive\Graduate study\Research\My_paper\Journal\Modal_ID_review\Code')
% Check data size
if size(y,1) < size(y,2)
    y = y';
end
if size(y_test,1) < size(y_test,2)
    y_test = y_test';
end

% Designate data size
[m, n] = size(y);
[m_test, n_test] = size(y_test);

if m ~= m_test
    TransformSize = m_test;
    fspan = linspace(0, fs, TransformSize);
else
    TransformSize = m;
    fspan = linspace(0, fs, TransformSize);
end

% FFT of the raw (or calibrated) modal test data
CH_all = 1:1:n;
CH_exclude = ~ismember(CH_all, CHref);
Yt = fft(y_test(:,CH_exclude),TransformSize);

% --------------------- Amplification via FRF ----------------------------
stds = std(y);    % Normalize the data (assuming Gaussian input)
means = mean(y);  % Obtain the mean of the data

% Data Normalization
yn = (y - means)./stds; % Normalize the data according to Gaussian first two moments

% Get the normalization factors of the reference channel
stds_ref = stds(CHref);
means_ref = means(CHref);

% Normalization factor for the response channels by getting rid of the
% reference channel normalization factors out of the other channels
means(CHref) = [];
stds(CHref) = [];

% Isolate the reference channel from the rest of the responses
ynref = yn(:,CHref);
yn(:, CHref) = [];

% ------------------ Data Inspection -------------------------------------
[Pyn, ~] = pwelch(yn,hanning(1024),512,[],fs);
[Pynref, freqs] = pwelch(ynref,hanning(1024),512,[],fs);

figure(1),clf
subplot(211)
for i = 1:size(yn,2)
pyn = plot(yn(:,i),'-k');
pyn.Color = [pyn.Color, 0.01*i];
hold on
end
pynref = plot(ynref,'r');
pynref.Color = [pynref.Color, 0.1];
title('Time Histories')
xlabel('Samples')
ylabel('Magnitude (V)')
subplot(212)
for i = 1:size(yn,2)
    p1 = plot(freqs, 10*log10(Pyn(:,i)), '--k');
    p1.Color = [p1.Color, 0.05*i];
    hold on
end
p2 = plot(freqs,10*log10(Pynref),'r');
p2.Color = [p2.Color, 1];
legend(p2,{'Reference'})
title('PSD Estimates')
xlabel('Frequency (Hz)')
ylabel('Magnitude/Frequency (dB/Hz)')
sgtitle('Calibration Data Set Visualization')

% ------------------------------------------------------------------------
%                   AMPLIFICATION SECTION BEGINS
% ------------------------------------------------------------------------
% FFT to the normalized data
Yn = fft(yn, TransformSize);
Ynref = fft(ynref, TransformSize);

% Channel-by-channel amplification
if isequal(lower(Method), 'tf')
    windsize = 256;
    noverlap = 256/2;
    % Define the frequency vector
    % Create an unevenly-spaced sample sequence
    fspan_base = linspace(0, fc, 2^10);
    fspan_frf = fspan_base.^2;
    indx = fspan_frf<fc;
    fspan_frf = fspan_frf(indx);
    fspan_frf = [fspan_frf fc];
    fspan = linspace(0, fc, TransformSize);
    % FRF estimation via Welch's method (FRF estimation takes time!!)
    [FRF, ~] = tfestimate(yn(:,1),yn(:,CHref),boxcar(windsize),noverlap,fspan_frf,fs,'Estimator','h1');
    FRFq = interp1(fspan_frf, FRF, fspan,'spline');
    clf
    plot(fspan_frf, mag2db(abs(FRF)), '*')
    hold on
    plot(fspan, mag2db(abs(FRFq)),'.')
    % Visualize the Amplifier
    figure
    % Concatenate the FRF to match the signal size
    FRFts = [FRFq, fliplr(FRFq)]';
    plot(abs(FRFts))
    % Amplify the original data using the FRF amplifier
    Yamp_FRF = Yn(:,1).*abs(FRFts);
    % Inverse FFT of the amplified data
    yamp_FRF = ifft(Yamp_FRF);
    yamp = real(yamp_FRF);
    yamp = stds(1)*(yamp + means(1));
elseif isequal(lower(Method), 'fft')
    % Define the frequency vector  
    FRF_FFT = zeros(size(Yn));         % Allocate memory for the FRF vector
    FRF = zeros(size(Yn));             % 
    Yamp_FRF = zeros(size(Yn));
    yamp_FRF = zeros(size(Yn));
    yamp = zeros(size(Yn));
    % This loop is for obtaining the FRF from the calibration data
    for i = 1:size(Yn,2)               % Loop through channels
        % FRF estimation via FFT
        FRF_FFT(:,i) = abs(Ynref)./abs(Yn(:,i)); % Assign FRF amplification factors to FRF_FFT
        FRF(:,i) = FRF_FFT(:,i);                 % Assign it to a new matrix
        % Detect inf/NaN case
        infindx = isinf(FRF(:,i));               % If zero freqency meets inf error
        FRF(infindx, i) = 0;                     % Adjust it to zero
        Yamp_FRF(:,i) = Yn(:,i).*FRF(:,i);       % Amplify the FFT(Y) according to FRF
        yamp_FRF(:,i) = ifft(Yamp_FRF(:,i));     % IFFT(Yamp_FRF), back to time domain
        yamp(:,i) = (yamp_FRF(:,i) + means(i))*stds(i); % Renormalizing the data
    end
    % This is for amplifying the modal test data y_test
    Yamp_test = zeros(size(Yt));
    yamp_test = zeros(size(Yt));
    for i = 1:size(Yt,2)               % Loop through channels
        Yamp_test(:,i) = Yt(:,i).*FRF(:,i);       % Amplify the test data according to FRF obtained from calibration data
        yamp_test(:,i) = ifft(Yamp_test(:,i));      % IFFT(Yamp_FRF), back to time domain
    end
end

% ------------------------------------------------------------------------
%                    OUTPUT VISUALIZATION SECTION BEGINS
% ------------------------------------------------------------------------

% Visualization to the FRF estiamte
figure(2),clf
plot(fspan, mag2db(abs(FRF)));
xlim([0 fc])
xlabel('Frequency (Hz)')
ylabel('FRF - $|H(f)|^{-1}$','Interpreter','latex')
title('Inverse FRF from FFT Estimates')
pbaspect([2 1 1])
% plot(fspan, FRF)

% Visualization of the amplification results
[Pyamp, ~] = pwelch(yamp,hanning(1024),512,[],fs);
figure(3),clf
subplot(211)
for i = 1:size(yn,2)
pyn = plot(yn(:,i),'-k');
pyn.Color = [pyn.Color, 0.01*i];
hold on
pyamp = plot(yamp(:,i),'-r');
pyamp.Color = [pyamp.Color, 0.01*i];
end
title('Time Histories')
xlabel('Samples')
ylabel('Magnitude (V)')
subplot(212)
for i = 1:size(yn,2)
    p1 = plot(freqs, 10*log10(Pyn(:,i)), '--k');
    p1.Color = [p1.Color, 0.05*i];
    hold on
    p2 = plot(freqs,10*log10(Pyamp(:,i)),'r');
    p2.Color = [p2.Color, 0.05*i];
end
legend([p1,p2],{'Original','Amplified'})
title('PSD Estimates')
xlabel('Frequency (Hz)')
ylabel('Magnitude/Frequency (dB/Hz)')
sgtitle('Amplified Calibration Results')

% Visualization of the amplification to the modal test data
[Pytest, ~] = pwelch(y_test,hanning(1024),512,[],fs);
[Pyatest, ~] = pwelch(yamp_test,hanning(1024),512,[],fs);
figure(4),clf
subplot(211)
for i = 1:size(yamp_test,2)
pyn = plot(y_test(:,i),'-k');
pyn.Color = [pyn.Color, 0.01*i];
hold on
pyamp = plot(yamp_test(:,i),'-r');
pyamp.Color = [pyamp.Color, 0.01*i];
end
title('Time Histories')
xlabel('Samples')
ylabel('Magnitude (V)')
subplot(212)
for i = 1:size(yn,2)
    p1 = plot(freqs, 10*log10(Pytest(:,i)), '--k');
    p1.Color = [p1.Color, 0.05*i];
    hold on
    p2 = plot(freqs,10*log10(Pyatest(:,i)),'r');
    p2.Color = [p2.Color, 0.05*i];
end
legend([p1,p2],{'Original','Amplified'})
title('PSD Estimates')
xlabel('Frequency (Hz)')
ylabel('Magnitude/Frequency (dB/Hz)')
sgtitle('Amplified Modal Test Data')


% ------------------------------------------------------------------------
%                              LEGACY CODE
% ------------------------------------------------------------------------

% --------------------- TEST INPUTS FOR DEBUGGING ------------------------
% Load reference data
% y = load('L:\My Drive\Graduate study\Research\My_paper\Journal\Modal_ID_review\Data\Base_acc_calibration_new.mat');
% y = y.data;
% if fc > 500
%     error('The maximum frequency range is 500Hz!')
% end
% windsize = 2^10;
% noverlap = windsize/2;
% Designate data size
% N = size(y,1);
% % FFT of the raw data
% Y = fft(y(:,2:end));
% ------------------------------------------------------------------------

%--------- Create the amplifier using curve fitting ---------
% a = 17.39;
% b = 0.01191;
% c = 62.78;
% % Frequency Variable
% x = 0:1:500;
% % Dynamic behavior function
% gf = b./(a*(exp(-(x-c)/a)+1));
% amp = max(gf) - gf;
% amp = gf/max(gf);
% amp = max(gf)./gf;
% % Describe the amplifier using piecewise polynomial
% amp_pp = spline(x,amp);
% 
% % plot the amplifier in frequency domain
% ff = linspace(0,fc,fix(size(data,1)/2));
% plot(ff, ppval(amp_pp, ff))
% FitAmp = ppval(amp_pp, ff);
% FitAmp = [FitAmp fliplr(FitAmp)]';
% Yamp_Fit = FitAmp.*Y(:,1);
% 
% % IFFT the amplified signal
% yamp_Fit = ifft(Yamp_Fit);

% FRF_spline = spline(f, FRF);
% %
% ff = linspace(0,fc,fix(size(data,1)/2));
% FitAmp_FRF = ppval(FRF_spline, ff);
% FitAmp_FRF = [FitAmp_FRF fliplr(FitAmp_FRF)]';