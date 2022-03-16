function [data_cal, caldata_cal, rel_mags] = SensorCalibrator(fs, NCH_cal, NCH_dat, RefIndx)
% Import the data
% ---------- Input parameters ----------
% Input: fs - Sampling frequency
%        NCH_cal - indices of channels used in the calibration data set
%        NCH_dat - indices of channels used in the modal test data set
%        RefIndx - index of the reference channel used for calibration
%        
% Output:
%        data_cal - calibrated testing data set
%        caldata_cal - calibrated calibration data set
%        rel_mags - calibration factor, relative standard deviations of
%        each channel signal with respect to the RefIndx channel.
%
% Hewenxuan Li, Oct 3, 2021
% ------------------------------------------------------------------------
% Use:
% This program can be used to calibrate the sensors that are subjected to
% variance biases between different sensors. First, a calibration data set
% using the collocation method is required. Second, a test data set is
% also needed and it will be calibrated by the calibration data set.
% Sampling frequency is used as a guide to plot the correct bandwidth of
% the selected data sets. The number of channels are desiganted to let the
% program know to which channels the calibration is conducted. Note that
% the calibration channels and the test channels HAVE TO BE IN THE SAME
% ORDER!
%
% By running this function, two open file windows will pop up in a roll.
% The first one asks for the CALIBRATION DATA SET, and the second one asks
% for the MODAL TEST DATA SET. After selecting the desired files, graphs
% that indicate the quality of the calibration will pop up.
%
% Example:
% Sampling frequency
% fs = 2000
% Number of channels in use
% NCH_cal = 6;                        % Number of active channels for calib
% NCH_dat = 6;                        
% RefIndx = 1;                        % Index of the reference channel
% ------------------------------------------------------------------------
% Define plot legend names
ChNames = {'CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8','CH9','CH9'};
warning("off")
% ----------------- Data Preprocessing -----------------------------------
fprintf('\n\n\n\n\n \n')
disp('============== Sensor Calibration Begins! ================')
% Add working path to openfile.m function
addpath("L:\My Drive\Graduate study\Research\Projects\Styles")
% Select the calibration data set
disp('Select the Calibration Data Set (from collocation method):')
[~, ~, CalibDataName] = openfile('*.*','- select the calibration data:');
caldata = readtable(CalibDataName{1});
disp(['Selected Calibration Data:', CalibDataName{1}])
% Select the modal testing data set
disp('Select the Data Set Need to be Adjusted:')
[~, ~, DataName] = openfile('*.*','- select the modal testing data:');
data = readtable(DataName{1});
disp(['Selected Modal Test Data:', DataName{1}])

% Convert to output type
caldata = table2array(caldata);
data = table2array(data);

% Check if the first channel is time channel
% Calibration data set
if mean(caldata(2:11,1) - caldata(1:10,1)) == 1/fs
    time = caldata(:,1);
    caldata(:,1) = [];
    disp('Time vector in calibration data detected and removed for amplification!')
    time_cal_flag = 1;
end
% Modal test data set
if mean(data(2:11,1) - data(1:10,1)) == 1/fs
    time_test = data(:,1);
    data(:,1) = [];
    disp('Time vector in modal test data detected and removed for amplification!')
    time_flag = 1;
end

data_original = data;
caldata_original = caldata;
[m, n] = size(data_original);
[mcal, ncal] = size(caldata_original);
% Trim data
caldata = caldata(:,NCH_cal);
data = data(:,NCH_dat);

% -------------------- Plot the raw data ---------------------------------
[pxx, f] = pwelch(caldata, hanning(1024), 512, [], fs);
figure(1),clf
subplot(221)
p1 = plot(f, 10*log10(pxx));
xlabel('Frequency (Hz)')
ylabel('Magnitude/Freq (dB/Hz)')
title('Power Spectral Density Estimate (Raw Calib. Data)')
pbaspect([2 1 1])
subplot(222)
for i = 1:size(caldata,2)
    [ksd, xi] = ksdensity(caldata(:,i),-5:0.1:5);
    histogram(caldata(:,i),-5:0.1:5,'Normalization','pdf','FaceAlpha',0.1,'EdgeAlpha',0.1)
    hold on
    plot(xi, ksd)
end
xlabel('Manitude (V)')
ylabel('Probability Density')
title('Probability Density Estimate (Raw Calib. Data)')
pbaspect([2 1 1])
legend(p1,{'CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8','CH9','CH9'})
hold off
% ------------ Calibration factors (means and stds) --------------------
% Select the first channel as the reference
caldata_n = zeros(size(caldata));          % noramlized data matrix initialization
for i = 1:size(caldata,2)          
    means(i) = mean(caldata(:,i));      % Mean subtraction
    caldata_n(:,i) = caldata(:,i) - mean(caldata(:,i));
    stds(i) = std(caldata(:,i));        % Normalize the mean-subtracted matrix
    caldata_n(:,i) = caldata(:,i)/stds(i); 
end
% Calculate the normalization factor as the ratio between the standard
% variance
% Relative magnitudes: = std(CH_i)/std(Ch_ref)
rel_mags = zeros(size(caldata,2),1); 
for i = 1:size(caldata,2)
    rel_mags(i) = stds(i)/stds(RefIndx); 
end
% Normalize (calibrate the sensors according to the variance/std)
caldata_cal = zeros(size(caldata));
for i = 1:size(caldata,2)
    caldata_cal(:,i) = ((caldata(:,i) - means(i))/rel_mags(i)) + means(i); % Normalization
end
% ------------------------- Plot the calibrated data ---------------------
[pxx, f] = pwelch(caldata_cal, hanning(1024), 512, [], fs);
subplot(223)
p1 = plot(f, 10*log10(pxx));
xlabel('Frequency (Hz)')
ylabel('Magnitude/Freq (dB/Hz)')
title('Power Spectral Density Estimate (Calibrated)')
pbaspect([2 1 1])
subplot(224)
for i = 1:size(caldata,2)
    [ksd, xi] = ksdensity(caldata_cal(:,i),-5:0.1:5);
    histogram(caldata_cal(:,i),-5:0.1:5,'Normalization','pdf','FaceAlpha',0.1,'EdgeAlpha',0.1)
    hold on
    plot(xi, ksd)
end
pbaspect([2 1 1])
xlabel('Manitude (V)')
ylabel('Probability Density Estimate')
title('Probability Density Estimate (Calibrated)')
legend(p1,{'CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8','CH9','CH9'})
hold off
sgtitle('Calibration Dataset Quality')
% ----------- Apply Calibration Information to the Raw Data -------------
data_cal = zeros(size(data));          % Allocate memory for the test data
for i = 1:size(caldata,2)          
    means_test(i) = mean(data(:,i));      % Mean subtraction
    data_cal(:,i) = data(:,i) - mean(data(:,i));
%     stds_test(i) = std(data(:,i));        % Normalize the mean-subtracted matrix
    data_cal(:,i) = data(:,i)/rel_mags(i) + means_test(i); % add mean value back
end

figure(2),clf
% ------------------------- Plot the calibrated data ---------------------
[pxx, f] = pwelch(data_cal, hanning(1024), 512, [], fs);
[pdata, ~] = pwelch(data, hanning(1024), 512, [], fs);
subplot(221)
p1 = plot(f, 10*log10(pdata));
pbaspect([2 1 1])
xlabel('Frequency (Hz)')
ylabel('Magnitude/Freq (dB/Hz)')
title('Power Spectral Density Estimate (raw test)')
subplot(222)
for i = 1:size(data,2)
    [ksd, xi] = ksdensity(data(:,i),-5:0.1:5);
    histogram(data(:,i),-5:0.1:5,'Normalization','pdf','FaceAlpha',0.1,'EdgeAlpha',0.1)
    hold on
    plot(xi, ksd)
end
pbaspect([2 1 1])
xlabel('Manitude (V)')
ylabel('Probability Density Estimate')
title('Probability Density Estimate (raw test)')
legend(p1,{'CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8','CH9','CH9'})
hold off

subplot(223)
p1 = plot(f, 10*log10(pxx));
pbaspect([2 1 1])
xlabel('Frequency (Hz)')
ylabel('Magnitude/Freq (dB/Hz)')
title('Power Spectral Density Estimate (Calibrated)')
subplot(224)
for i = 1:size(data_cal,2)
    [ksd, xi] = ksdensity(data_cal(:,i),-5:0.1:5);
    histogram(data_cal(:,i),-5:0.1:5,'Normalization','pdf','FaceAlpha',0.1,'EdgeAlpha',0.1)
    hold on
    plot(xi, ksd)
end
pbaspect([2 1 1])
xlabel('Manitude (V)')
ylabel('Probability Density Estimate')
title('Probability Density Estimate (Calibrated)')
legend(p1,{'CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8','CH9','CH9'})
hold off
sgtitle('Calibrated Modal Test Data')

% ----------------------- Output data curing -----------------------------
if time_cal_flag == 1 && time_flag == 1
    time_test = time_test';                    % Time vector for test data
    time_cal = time';                     % Time vector for calibration data
else
    time_test = 0:1/fs:1/fs*(m - 1);          % Time vector for test data
    time_cal = 0:1/fs:1/fs*(mcal - 1);   % Time vector for calibration data
end
CH_all = 1:1:n;                                % All channel indices
NCH_excluded = ~ismember(CH_all,NCH_dat);      % Channels exluded from the calibration
NCH_excluded(1) = false;                       % Negate the first uncalibrated channel (time vector)
data_cal = [time_test', data_cal, data_original(:,NCH_excluded)]; % Padding time vector and unchanged channel signals
caldata_cal = [time_cal', caldata_cal, caldata_original(:,NCH_excluded)]; % Padding time vector and unchanged channel signals

disp('============== Sensor Calibration Done! ================')
return