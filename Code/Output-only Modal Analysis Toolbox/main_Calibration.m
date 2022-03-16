% ========================== Calibration =================================
% ------------------------------------------------------------------------
%                      DATA CALIBRATION SECTION
% ------------------------------------------------------------------------
%
%
% As part of the demo for the output-only modal analysis toolbox
% Work in progress...
% Hewenxuan Li, Oct 2021
%% Without data cleasing
% ------------------------------------------------------------------------
%              DATA CALIBRATION SECTION - NO AMPLIFICATION
% ------------------------------------------------------------------------
% DAQ information
fs = 2000;      % Sampling frequency
fc = 1000;       % Cutoff frequency
NCH_cal = 1:8;  % Channel number for the calibration data set
NCH_dat = 1:8;  % Channel number for the modal test data set
RefIndx = 1;    % Reference channel used in the calibration data set

% ----------------------------- EXMALE DATA ------------------------------
% Choose the following data sets as illustration:
% 1. Calibration_333B42_WGN_fs_1000_fc_500_500mvpp_002.txt
% 2. Calibration_333B42_Harmonic4_fs_1000_fc_500_500mvpp_001.txt
% ------------------------------------------------------------------------

[data_cal, caldata_cal, rel_mags] = SensorCalibrator(fs, NCH_cal, NCH_dat, RefIndx);

% ------------------------------------------------------------------------
%                      FRF AMPLIFICATION SECTION
% ------------------------------------------------------------------------

% Method = 'FFT';
% y = caldata_cal(:,2:end);
% y_test = data_cal(:,2:end);
% CHref = 8;
% 
% [yamp_test, FRF] = AccAmplifier(y, y_test, CHref, fs, fc, Method);
% Assign the calibrated data to a new variable

y = caldata_cal;
y_test = data_cal;
CHref = 8;         % REFERENCE CHANNEL
CHamp = 1:6;       % RESPONSE CHANNELS
% Run amplification
[yamp_test, FRF] = AccAmplifier_TF(y, y_test, CHref, CHamp, fs, fc, 'linear', 2^10, 1024, [1 400], 'Preamp');

% ------------------------------------------------------------------------

%% With data cleansing
% ------------------------------------------------------------------------
%              DATA CALIBRATION SECTION - WITH AMPLIFICATION
% ------------------------------------------------------------------------
% DAQ information
addpath("L:\My Drive\Graduate study\Research\Projects\OS")
fs = 2000;      % Sampling frequency
fc = 1007;       % Cutoff frequency
NCH_cal = 1:7;  % Channel number for the calibration data set
NCH_dat = 1:7;  % Channel number for the modal test data set
RefIndx = 1;    % Reference channel used in the calibration data set

% ----------------------------- EXMALE DATA ------------------------------
% Choose the following data sets as illustration:
% 1. Calibration_333B42_WGN_fs_1000_fc_500_500mvpp_002.txt
% 2. Calibration_333B42_Harmonic4_fs_1000_fc_500_500mvpp_001.txt
% ------------------------------------------------------------------------

[data_cal, caldata_cal, rel_mags] = SensorCalibrator(fs, NCH_cal, NCH_dat, RefIndx);

% ------------------------------------------------------------------------
%                     FRF AMPLIFICATION SECTION
% ------------------------------------------------------------------------

% Method = 'FFT';
% y = caldata_cal(:,2:end);
% y_test = data_cal(:,2:end);
% CHref = 8;
% 
% [yamp_test, FRF] = AccAmplifier(y, y_test, CHref, fs, fc, Method);

y = caldata_cal;
y_test = data_cal;
CHref = 8;         % REFERENCE CHANNEL
CHamp = 1:6;       % RESPONSE CHANNELS
InterpMethod = 'linear'; % Method used for interpolating the frf
Nfreq = 2^11;      % # of points used for samping the frf
WindSize = 2^10;   % Window size for frf estimation
DeadBand = [0.1 800]; % Dead band boundaries
AmpMethod = 'Preamp'; % Amplification strategy
% Run amplification (calling AccAmplifier_TF.m)
[yamp_test, FRF] = AccAmplifier_TF(y, y_test, CHref, CHamp, fs, fc, InterpMethod, Nfreq, WindSize, DeadBand, AmpMethod);
data = [];
data = [y_test(:,1), yamp_test, y_test(:,8:9)];

% ------------------------------------------------------------------------
%                LOADING DATA SUBROUTINE FOR SAVING INFO
% ------------------------------------------------------------------------      
% Basic inforamtion 
InputIndx = 7; % 
ExpWindRate = 0.7;
ExpWindSize = 0.7;
% Select the data set
% ------------------------------------------------------------------
% Choose the following files for illustration:
% ..\Impulse_333B42_fs_1000_fc_500_N_5_i_3_001.txt
% ------------------------------------------------------------------
addpath("L:\My Drive\Graduate study\Research\Projects\OS")
disp('Select the Data Set Need to be Adjusted:')
[file, path, DataName] = openfile('*.*','- select the modal testing data again for saving the data:');
SaveName = split(file{1},'.');
SaveName = SaveName{1};

DataInspection(data)

% DATA CLEANSING SUBROUTINE
[data_amped_cln_out, windsize] = DataCleansing(data, fs, InputIndx, NCH_dat, ExpWindRate, ExpWindSize);

DataInspection(data_amped_cln_out)

impulse_amped_cleansed_data.data = data_amped_cln_out;
impulse_amped_cleansed_data.windsize = windsize;
impulse_amped_cleansed_data.fs = fs;
save([path, SaveName, '_amped_cln.mat'], 'impulse_amped_cleansed_data')