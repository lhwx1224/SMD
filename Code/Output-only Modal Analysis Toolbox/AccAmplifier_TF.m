function [yamp_test, FRFq] = AccAmplifier_TF(y, y_test, CHref, CHamp, fs,...
    fc, FRF_scale, Nf, windsize, deadband, DeadBandMode)
% AccAmplifier_TF amplifies the modal test data 'y_test' based on a set of
% calibration data set 'y' (which may be calibrated first). A reference
% channel 'CHref' is used to designate the targeted channel signal based 
% off of which the amplifier is designed. Sampling parameters, sampling
% frequency 'fs' and cut-off frequency 'fc', are required for generating
% plots. 'Method' is used to choose from either Welch's method or FFT
% method for estimating the frequency response (FRF) functions.
% 
% Syntax: [yamp, yamp_FRF] = AccAmplifier(y, y_test, CHref, fs, fc,
% FRF_scale, Nf, windsize, deadband);
% 
% Input: y - 2d array, calibration data set
%        y_test - 2d array, modal test data set
%        CHref - int, reference channel index
%        CHamp - int or array, index of channels that are being amplified
%                IF int, put CHamp = 1, this impose only ONE FRF amplifier
%                for ALL channels;
%                IF an array, e.g., CHamp = [1, 3, 7], it only amplifies
%                channel index 1, 3, and 7 for the given data set based on
%                their corresponding channels used in the calibration data
%                set.
%        fs - double, sampling frequency
%        fc - double, cut-off frequency
%        FRF_scale - string, frequency vector used to sample the frequency
%                    response function estimate
%                    (1) 'Linear' - linear f vector for FRF
%                    (2) 'Nonlinear' - nonlinear f vector for FRF
%        Nf - number of frequency points used for generating the FRF scale
%        windsize - int, window size for the FRF estimate
%        deadband - 2-by-1 array, band-pass range indicate the deadband
%
% Output: yamp_test - 2d array, amplified modal test data
%         FRFq - 2d array, frequency response function estimate for
%                    amplification 
%
% Required functions: DataInspection.m 
% Hewenxuan Li, Created on April 23, 2021
% Last Modified: 
% October 1, 2021: Added FFT-based amplitude amplification
% October 4, 2021: fully functional to directly incorporate both
% calibration data set and modal testing data set.
% October 6, 2021: based on the AccAmplifier.m function, this code is now
% dedicated to use the transfer function estimate. New variables are added
% to control (1) the frequency variable to control the sample point
% density, (2) the number of frequency sample points used in the FRF/TF
% estimation, (3) the window size used for FRF/TF estimate, and (4) the
% deadband of the FRF/TF estimate to avoid excessive undulations during the
% estimation.
% October 25, 2021: added more option to the CHamp, now it allows all
% channel data use only one specific FRF based on CHamp and CHref. For
% example, when CHamp = 2 and CHref = 1, the FRF amplifier is based on the
% Channel 2 over Channel 1 signals. 
disp('=================== AccAmplifier_TF() Begins: =====================')
% Add toolbox path
addpath('L:\My Drive\Graduate study\Research\My_paper\Journal\Modal_ID_review\Code')
% ----------------------- CHECK DATA SIZE --------------------------------
if size(y,1) < size(y,2)
    y = y'; % Calibration data
end
if size(y_test,1) < size(y_test,2)
    y_test = y_test'; % Modal test data
end

% ------------------------ IDENTIFY TIME VECTORS -------------------------
if mean(y(2:11,1) - y(1:10,1)) == 1/fs
    y_time = y(:,1);      % Calibration data set
    y(:,1) = [];
    disp('Time vector in calibration data detected and removed for amplification!')
end

if mean(y_test(2:11,1) - y_test(1:10,1)) == 1/fs
    y_time_test = y_test(:,1);  % Modal test data set
    y_test(:,1) = [];
    disp('Time vector in modal test data detected and removed for amplification!')
end

% ----------------- GET DATA SET SIZES -----------------------------------
[m, n] = size(y);
[m_test, n_test] = size(y_test);
% Change data sizes to prepare for the amplification
if mod(m_test,2) ~= 0
    y_test = y_test(1:end-1,:);
    [m_test, n_test] = size(y_test);
end

if mod(m,2) ~= 0
    y = y(1:end-1,:);
    [m, n] = size(y);
end
% Assign size of the Welch's FRF estimates
if m ~= m_test % IF the test data set has different size from the calibration set
    TransformSize = m_test;                       % Identify the FRF size (same size as the temporal stamps)
    fspan = linspace(0, fc, TransformSize/2);     % Frequency vector (test data) designated by the cutoff frequency
    fspanc = linspace(0, fc, m/2);                % Frequency vector (calibration data) designated by ...
else % IF They have the same sizes
    TransformSize = m;                                            
    fspan = linspace(0, fc, TransformSize/2);
end

% FFT of the raw (or calibrated) modal test data (exclude the reference
% channel from the rest of the response data) 
% CHall = 1:1:n;
% CHamp = 1:6;   % Designate the channel indices to be amplified
% CH_exclude = ~ismember(CHall, CHamp); % excluded channels for amplification
% Yt = fft(y_test(:,CHamp),TransformSize);   % FFT(y_test) frequency spectrum of the test data

% --------------------- Calibration Data Setting -------------------------
stds = std(y);    % Normalize the data (assuming Gaussian input)
means = mean(y);  % Obtain the mean of the data

% Normalized and time-vector-excluded data set
yn = (y - means)./stds; % Normalize the data according to Gaussian first two moments

% Isolate the reference channel from the rest of the responses
ynref = yn(:,CHref);             % Normalized reference data
yn(:, CHref) = [];               % Noramlized calibration data

% --------------- Normalization factors from calibration data ------------
% Get the normalization factors of the reference channel
stds_ref = stds(CHref);
means_ref = means(CHref);
% Normalization factor for the response channels by getting rid of the
% reference channel normalization factors out of the other channels
means(CHref) = [];
stds(CHref) = [];

% ------------------ FFT to all the prepared data sets -------------------
if length(CHamp) == 1                % IF CHamp is an upper bound
    Yt = fft(y_test(:,:));     % FFT to the modal test data based on all channels
else                                 % ELSEIF CHamp is a vector
    Yt = fft(y_test(:,CHamp));       % FFT to the modal test data
end
Yn = fft(yn);                    % FFT to the normalized calibration data
Ynref = fft(ynref);              % FFT to the normalized reference calibration data (NOT IN USE, CONSIDER ERASING)

% ================== Data Inspection =====================================
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
% ========================================================================

% ------------------------------------------------------------------------
%                   AMPLIFICATION SECTION SETUP
% ------------------------------------------------------------------------
% FTF Estimation Parameters
if nargin < 8
    Nf = 2^10;                   % number of frequency sample points
    windsize = 1024;             % Window Size (default 1024)
    deadband = [];
    DeadBandMode = 'postamp';
elseif nargin < 9
    windsize = 1024;
    deadband = [];
    DeadBandMode = 'postamp';
elseif nargin < 10
    deadband = [];
    DeadBandMode = 'postamp';
elseif nargin < 11
    DeadBandMode = 'postamp';
end
noverlap = windsize/2;           % Number of points used for overlapping

% ---------------- Define the frequency vector ---------------------------
% Create an unevenly-spaced sample sequence, there are only Nf query points
% used for the esimtation to (1) accelerate the calculation and (2) to
% promote the averaging... 

% IF using the nonlinear FRF points
if isequal(lower(FRF_scale), 'nonlinear') % IF NONLINEAR FREQ SCALE USED
    fspan_base = linspace(0, fc, Nf); % Create linearly-spaced base frequency
    fspan_frf = fspan_base.^3;        % Create nonlinearly-spaced frequency points
    indx = fspan_frf<fc;              % If the range is less than cutoff,
    fspan_frf = fspan_frf(indx);      % padding the fc to the end of the freqs
    fspan_frf = [fspan_frf fc]';
elseif isequal(lower(FRF_scale), 'linear') % IF LINEAR FREQ SCALE USED
    % IF using the linear FRF points
    fspan_frf = linspace(0, fc, Nf);  % Linear range
end
% Allocate memory for the estimation
% FRF estimation via Welch's method (FRF estimation takes most of the time!!)
FRF = zeros(length(fspan_frf),length(CHamp));  % Nf*# of CHamp freq. vector (for FRF estimation)  
FRFq = zeros(TransformSize/2, length(CHamp));  % # of Half freq*# of CHamp new query points for full-size transformation
FRFts = zeros(TransformSize, length(CHamp));   % # of Full freq*# of CHamp new query points for IFFT

% ------------------------------------------------------------------------
%                     AMPLIFICATION SECTION BEGINS
% Two modes can be selected:
% 1. pre-amplification - amplifies the signal and apply the deadband
% adjustment before interpolation.
% 2. post-amplification - amplifies the signal and interpolate the FRF
% first; then, conduct deadband adjustment afterwards.
% ------------------------------------------------------------------------
if isequal(lower(DeadBandMode), 'preamp') % IF preamplification is used
disp('Pre-amplification is in use! FRF estimation Begins!')
    for i = 1:length(CHamp)
        [frf, frf_f] = tfestimate(yn(:,CHamp(i)), ynref ,boxcar(windsize), noverlap, fspan_frf, fs,'Estimator','h1');
        % --------- NEW TESTING CODE 10/7/21 ---------------------------
        if ~isempty(deadband)
            disp(['Deadband is applied from 0-',num2str(deadband(1)),' and from ',num2str(deadband(2)),'-',num2str(fc),'.'])
            indxfrf = frf_f > deadband(2) | frf_f < deadband(1);
            FRF(:,i) = frf;
            FRF(indxfrf,i) = 1;
        else
            disp('Full FRF estimate is used without deadbands!')
            FRF(:,i) = frf;
        end
        % --------------------------------------------------------------
        FRFq(:,i) = interp1(fspan_frf, FRF(:,i), fspan,'linear')';
        FRFqc(:,i) = interp1(fspan_frf, FRF(:,i), fspanc,'linear')';
        FRFts(:,i) = [FRFq(:,i); flipud(FRFq(:,i))];
        FRFtsc(:,i) = [FRFqc(:,i); flipud(FRFqc(:,i))];
    end

else % Elseif post amplification with deadband
    disp('Post-amplification is in use! FRF estimation Begins!')
    for i = 1:length(CHamp)
        [frf, frf_f] = tfestimate(yn(:,CHamp(i)), ynref ,boxcar(windsize), noverlap, fspan_frf, fs,'Estimator','h1');
        FRF(:,i) = frf;
        % --------------------------------------------------------------
        FRFq(:,i) = interp1(fspan_frf, FRF(:,i), fspan,'spline')';
        FRFqc(:,i) = interp1(fspan_frf, FRF(:,i), fspanc,'spline')';
        FRFts(:,i) = [FRFq(:,i); flipud(FRFq(:,i))];
        FRFtsc(:,i) = [FRFqc(:,i); flipud(FRFqc(:,i))];
    end
    % ---------------- 10/6/2021 MOD ---------------------------
    if ~isempty(deadband)
        disp(['Deadband is applied from 0-',numstr(deadband(1)),' and from ',num2str(deadband(2)),'-',num2str(fc),'.'])
        % Deadband assignment
        indxq = fspan > deadband(2) | fspan < deadband(1);
        indxqc = fspanc > deadband(2) | fspanc < deadband(1);
        for i = 1:length(CHamp)
            FRFq(:,i) = abs(FRFq(:,i));
            FRFqc(:,i) = abs(FRFqc(:,i));
            FRFq(indxq,i) = 1;
            FRFqc(indxqc,i) = 1;
            FRFts(:,i) = [FRFq(:,i); flipud(FRFq(:,i))];
            FRFtsc(:,i) = [FRFqc(:,i); flipud(FRFqc(:,i))];
        end
    else
        disp('Full FRF estimate is used without deadbands!')
    end
% ---------------- END MOD ---------------------------------
end % END if condition on the mode of applying the deadband!

% ------------------------------------------------------------------------
%                    Channel-by-channel amplification
%  IF only one channel of calibration data is used for amplification, all
%  the modal test data will be amplified using the same FRF estimate;
%  ELSE each channel will be amplified separately based on their
%  corresponding calibration data counterparts.
% ------------------------------------------------------------------------
disp(['Applying the amplifier to the following channels:', num2str(CHamp), '.'])
% Amplify the original data using the FRF amplifier
Yamp_test = zeros(size(Yt));
yamp_test = zeros(size(Yt));

Ynamp = zeros(size(Yn));
ynamp = zeros(size(Yn));
% --------------- MOD 10/25/21 Begins here -----------
if length(CHamp) == 1
    for i = 1:size(Yt, 2)
        % --------- The calibration data set ------
        Ynamp(:,i) = Yn(:,1).*abs(FRFtsc(:,1));
        ynamp(:,i) = real(ifft(Ynamp(:,i)));
        % --------- The modal test data set -------
        Yamp_test(:,i) = Yt(:,i).*abs(FRFts(:,1));    % Amplify using only one channel
        yamp_test(:,i) = real(ifft(Yamp_test(:,i)));  % IFFT to obtain the amplified signal
    end
else
% --------------- MOD 10/25/21 Ends here -------------
    for i = 1:size(Yt, 2)
        % --------- The calibration data set ------
        Ynamp(:,i) = Yn(:,1).*abs(FRFtsc(:,i));
        ynamp(:,i) = real(ifft(Ynamp(:,i)));
        % --------- The modal test data set -------
        Yamp_test(:,i) = Yt(:,i).*abs(FRFts(:,i));    % Amplify using each channel respectively
        yamp_test(:,i) = real(ifft(Yamp_test(:,i)));  % IFFT to obtain the amplified signal
    end
end
% Inverse FFT of the amplified data

% ------------------------------------------------------------------------
%                    OUTPUT VISUALIZATION SECTION BEGINS
% ------------------------------------------------------------------------
disp('Amplification finished! Outputing visualization...')
% Visualization of the raw FRF and the interpolated FRF
figure(1),clf
plot(fspan_frf, mag2db(abs(FRF)), '*')
hold on
plot(fspan, mag2db(abs(FRFq)),'.')
pbaspect([2 1 1])
xlabel('Frequency (Hz)')
ylabel('FRF - $|H(f)|^{-1}$','Interpreter','latex')
title('Inverse FRF from TF Estimates')
% Visualization to the FRF estiamte


% Visualization of the amplification results
[Pynamp, ~] = pwelch(ynamp,hanning(1024),512,[],fs);
figure(3),clf
subplot(211)
if length(CHamp) == 1
    for i = 1:size(yn, 2)
        pyn = plot(yn(:,i),'-k');
        pyn.Color = [pyn.Color, 0.01*i];
        hold on
        pyamp = plot(ynamp(:,i),'-r');
        pyamp.Color = [pyamp.Color, 0.01*i];
    end
else
    for i = 1:length(CHamp)
        pyn = plot(yn(:,CHamp(i)),'-k');
        pyn.Color = [pyn.Color, 0.01*i];
        hold on
        pyamp = plot(ynamp(:,i),'-r');
        pyamp.Color = [pyamp.Color, 0.01*i];
    end
end
title('Time Histories')
xlabel('Samples')
ylabel('Magnitude (V)')
subplot(212)
if length(CHamp) == 1
    for i = 1:size(yn,2)
        p1 = plot(freqs, 10*log10(Pyn(:,i)), '--k');
        p1.Color = [p1.Color, 0.05*i];
        hold on
        p2 = plot(freqs,10*log10(Pynamp(:,i)),'r');
        p2.Color = [p2.Color, 0.05*i];
    end 
else
    for i = 1:length(CHamp)
        p1 = plot(freqs, 10*log10(Pyn(:,i)), '--k');
        p1.Color = [p1.Color, 0.05*i];
        hold on
        p2 = plot(freqs,10*log10(Pynamp(:,i)),'r');
        p2.Color = [p2.Color, 0.05*i];
    end
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
pyn.Color = [pyn.Color, 0.1*i];
hold on
pyamp = plot(yamp_test(:,i),'-r');
pyamp.Color = [pyamp.Color, 0.1*i];
end
title('Time Histories')
xlabel('Samples')
ylabel('Magnitude (V)')
subplot(212)
if length(CHamp) == 1
    for i = 1:size(Pyatest, 2)
        p1 = plot(freqs, 10*log10(Pytest(:,i)), '--k');
        p1.Color = [p1.Color, 0.1*i];
        hold on
        p2 = plot(freqs,10*log10(Pyatest(:,i)),'r');
        p2.Color = [p2.Color, 0.1*i];
    end
else
    for i = 1:length(CHamp)
        p1 = plot(freqs, 10*log10(Pytest(:,i)), '--k');
        p1.Color = [p1.Color, 0.1*i];
        hold on
        p2 = plot(freqs,10*log10(Pyatest(:,i)),'r');
        p2.Color = [p2.Color, 0.1*i];
    end
end
legend([p1,p2],{'Original','Amplified'})
title('PSD Estimates')
xlabel('Frequency (Hz)')
ylabel('Magnitude/Frequency (dB/Hz)')
sgtitle('Amplified Modal Test Data')

disp('=================== AccAmplifier_TF() Done! =======================')

return