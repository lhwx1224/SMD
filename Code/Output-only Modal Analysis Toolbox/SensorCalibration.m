%% Import the data
% ---------- Input parameters ----------
% Sampling frequency
fs = 2000;
% Number of channels in use
NCH_cal = 6;                        % Number of active channels for calib
NCH_dat = 6;                        
refindx = 1;                        % Index of the reference channel

% Define plot legend names
ChNames = {'CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8','CH9','CH9'};

% ----------------- Data Preprocessing -----------------------------------
% Add working path to openfile.m function
addpath("L:\My Drive\Graduate study\Research\Projects\Styles")
% Select the calibration data set
disp('Select the Calibration Data Set (from collocation method):')
[~, ~, CalibDataName] = openfile('*.*','- select the calibration data:');
caldata = readtable(CalibDataName{1});

% Select the modal testing data set
disp('Select the Data Set Need to be Adjusted:')
[~, ~, DataName] = openfile('*.*','- select the modal testing data:');
data = readtable(DataName{1});

% Convert to output type
caldata = table2array(caldata);
data = table2array(data);

% Trim data
caldata = caldata(:,2:NCH_cal + 1);
data = data(:,2:NCH_dat + 1);

% -------------------- Plot the raw data ---------------------------------
[pxx, f] = pwelch(caldata, hanning(1024), 512, [], fs);
figure(1),clf
subplot(221)
p1 = plot(f, 10*log10(pxx));
xlabel('Frequency (Hz)')
ylabel('Magnitude/Freq (dB/Hz)')
title('Power Spectral Density Estimate (Raw Calib. Data)')
subplot(222)
for i = 1:size(caldata,2)
    [f, xi] = ksdensity(caldata(:,i),-5:0.1:5);
    histogram(caldata(:,i),-5:0.1:5,'Normalization','pdf','FaceAlpha',0.1,'EdgeAlpha',0.1)
    hold on
    plot(xi, f)
end
xlabel('Manitude (V)')
ylabel('Probability Density')
title('Probability Density Estimate (Raw Calib. Data)')
legend(p1,{'CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8','CH9','CH9'})
hold off
%% Calibration factors (means and stds)
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
    rel_mags(i) = stds(i)/stds(refindx); 
end
% Normalize (calibrate the sensors according to the variance/std)
caldata_cal = zeros(size(caldata));
for i = 1:size(caldata,2)
    caldata_cal(:,i) = caldata(:,i)/rel_mags(i); % Normalization
end
% ------------------------- Plot the calibrated data ---------------------
[pxx, f] = pwelch(caldata_cal, hanning(1024), 512, [], 2000);
subplot(223)
p1 = plot(f, 10*log10(pxx));
xlabel('Frequency (Hz)')
ylabel('Magnitude/Freq (dB/Hz)')
title('Power Spectral Density Estimate (Calibrated)')
subplot(224)
for i = 1:size(caldata,2)
    [f, xi] = ksdensity(caldata_cal(:,i),-5:0.1:5);
    histogram(caldata_cal(:,i),-5:0.1:5,'Normalization','pdf','FaceAlpha',0.1,'EdgeAlpha',0.1)
    hold on
    plot(xi, f)
end
xlabel('Manitude (V)')
ylabel('Probability Density Estimate')
title('Probability Density Estimate (Calibrated)')
legend(p1,{'CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8','CH9','CH9'})
hold off
%% Apply Calibration Information to the Raw Data
data_n = zeros(size(data));          % Allocate memory for the test data
for i = 1:size(caldata,2)          
    means_test(i) = mean(data(:,i));      % Mean subtraction
    data_n(:,i) = data(:,i) - mean(data(:,i));
%     stds_test(i) = std(data(:,i));        % Normalize the mean-subtracted matrix
    data_n(:,i) = data(:,i)/rel_mags(i) + means_test(i); % add mean value back
end

figure(2),clf
% ------------------------- Plot the calibrated data ---------------------
[pxx, f] = pwelch(data_cal, hanning(1024), 512, [], 2000);
subplot(121)
p1 = plot(f, 10*log10(pxx));
xlabel('Frequency (Hz)')
ylabel('Magnitude/Freq (dB/Hz)')
title('Power Spectral Density Estimate (Calibrated)')
subplot(122)
for i = 1:size(data_cal,2)
    [f, xi] = ksdensity(data_cal(:,i),-5:0.1:5);
    histogram(data_cal(:,i),-5:0.1:5,'Normalization','pdf','FaceAlpha',0.1,'EdgeAlpha',0.1)
    hold on
    plot(xi, f)
end
xlabel('Manitude (V)')
ylabel('Probability Density Estimate')
title('Probability Density Estimate (Calibrated)')
legend(p1,{'CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8','CH9','CH9'})
hold off
