% ------------------------------------------------------------------------
%               MAIN FILE OF THE DATA CLEASING FUNCTION
% ------------------------------------------------------------------------
%
%
% As part of the demo program for output-only modal analysis toolbox
% Working in progress...
% 
% Hewenxuan Li, Oct 2021

%% LOADING DATA
% Set up the working environment 
addpath("L:\My Drive\Graduate study\Research\Projects\OS")
% Indices of channels used for calib and amp from the modal test
NCH_dat = [1:7];       
% Basic inforamtion 
fs = 2000;
InputIndx = 8; %
SaveType = '*.mat';
ExpWindRate = 0.5;
ExpWindSize = 0.5;
% Select the data set
% ------------------------------------------------------------------
% Choose the following files for illustration:
% ..\Impulse_333B42_fs_1000_fc_500_N_5_i_3_001.txt
% ------------------------------------------------------------------
disp('Select the Data Set Need to be Adjusted:')
[file, path, DataName] = openfile('*.*','- select the modal testing data:');
data = readtable(DataName{1});
data = table2array(data);

SaveName = split(file{1},'.');
SaveName = SaveName{1};

DataInspection(data)
%% DATA CLEANSING
[data_cln_out, windsize] = DataCleansing(data, fs, InputIndx, NCH_dat, ExpWindRate, ExpWindSize);

DataInspection(data_cln_out)

%% DATA SAVING
impulse_cleansed_data.data = data_cln_out;
impulse_cleansed_data.windsize = windsize;
impulse_cleansed_data.fs = fs;
save([path, SaveName, '_cln.mat'], 'impulse_cleansed_data')