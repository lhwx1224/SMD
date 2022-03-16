%% Data Pre-processing for the impact modal testing
% (1) Data calibration
% (2) Data Amplification (optional)
% (3) Data Cleansing (windowing and truncating)
% (4) Modal Identification for Pre-processing (if the data is good)
% (5) Data Saving
%% SELECT DATA SETS
addpath("L:\My Drive\Graduate study\Research\Projects\OS")
[file, path, selectedfile] = openfile('*.txt', 'Select the modal data file:');
FileName = split(file,'.');
FileName = FileName{1};
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import and concatenate the data
data = cell(length(selectedfile),1);
for i = 1:length(selectedfile)
    data{i} = table2array(readtable(selectedfile{i}, opts));
end

data = cell2mat(data);
%% DATA PREPROCESSING
% ------------------------------------------------------------------------
%                PART I - SENSOR CALIBRATION
% ------------------------------------------------------------------------
% Setup Environment
addpath("L:\My Drive\Graduate study\Research\Projects\OS")
Amp = 0;        % Use amplification?
fs = 1000;      % Sampling frequency
fc = 500;      % Cutoff frequency
NCH_cal = 1:7;  % Channel number for the calibration data set
NCH_dat = 1:7;  % Channel number for the modal test data set
RefIndx = 1;    % Reference channel used in the calibration data set

[data_cal, caldata_cal, rel_mags] = SensorCalibrator(fs, NCH_cal, NCH_dat, RefIndx);

% ------------------------------------------------------------------------
%                PART II - SENSOR AMPLIFICATION
% ------------------------------------------------------------------------
if Amp
y = caldata_cal;
y_test = data;
CHref = 8;         % REFERENCE CHANNEL
CHamp = 1:7;       % RESPONSE CHANNELS
CHbyCHAmp = 1;
% [yamp_test, FRF] = AccAmplifier_TF(y, y_test, CHref, CHamp, fs, fc, 'linear', 2^10, 1024, [1 900], 'Preamp');
[yamp_test, FRF] = AccAmplifier_TF_Edit(y, y_test, CHref, CHamp, fs, fc, 'linear', 2^10, 1024, [1 900], 'Preamp', CHbyCHAmp);
% Pending Time vector to the first colulmn of the amplified data
yamp_test = [data(1:size(yamp_test,1),1), yamp_test, data(1:size(yamp_test,1),9)];
else
    yamp_test = data;
end
% ------------------------------------------------------------------------
%                PART III -  DATA CLEANSING
% ------------------------------------------------------------------------
InputIndx = 8;
NCH_dat = 1:6;
ExpWindRate = 0.2;
ExpWindSize = 1;
% Data cleansing
[data_cln_out, windsize] = DataCleansing(yamp_test, fs, InputIndx, NCH_dat, ExpWindRate, ExpWindSize, 'False', 50);

% Assign the pre-processed data
data_prep = data_cln_out;
% 
tfestimate(data_prep(:,1), data_prep(:,2:7), boxcar(windsize),0,[],2000,'Estimator','H1');
set(gcf, 'Papersize', [6 6*0.5])
set(gcf, 'PaperPosition', [0 0 6 6*0.5])

%% Modal FRF esitmates
% ------------------------------------------------------------------------
%             PART IV - MODAL IDENFICATION BEGINS (PRE)
% ------------------------------------------------------------------------
nm = 4;
Nsensors = 6;
[FRF, f] = modalfrf(data_prep(:,1), data_prep(:,2:7), fs, boxcar(windsize),0, 'Sensor', 'dis');

freq_range = {[0, 50],[50, 100],[100, 300],[400,500]};
ms = zeros(Nsensors,nm);
fn = zeros(nm,1);
dr = zeros(nm,1);
% fn - natural frequencies, dr - damping ratios, ms - mode shapes
for i = 1:length(freq_range)
    [fn_temp, dr_temp, ms_temp] = modalfit(FRF, f, fs, 1, 'FitMethod', 'lsrf', 'FreqRange', freq_range{i});
    fn(i) = fn_temp;
    dr(i) = dr_temp;
    ms(:,i) = ms_temp;
end

Phi = real(ms);
Phi_norm = Phi./Phi(1,:);
Phi_norm = [Phi_norm; zeros(1,size(Phi_norm,2))];
figure(1),clf
p1 = plot(Phi_norm);
hold on
plot([1 size(Phi_norm,1)],[0 0],'--','Color','#B6B6B2')
pbaspect([2 1 1])
title('Experimental Modal Analysis --- LSRF')
xlabel('Sensor Location')
ylabel('Normalized Magnitude')
grid on
% legend(p1,ModeLegends,'Location','EastOutside','NumColumns',2)
% title('Cleansed-and-amplified Data Set')
set(gcf, 'Papersize', [6 6*0.5])
set(gcf, 'PaperPosition', [0 0 6 6*0.5])


%% SAVE DATA TO MAT FORMAT
% ------------------------------------------------------------------------
%            PART V - SAVE DATA TO .MAT FORMAT
% ------------------------------------------------------------------------
SaveName = [path, FileName, '.mat'];
save(SaveName, 'data_cln_out')
disp('------------------ Mat file saved ! ---------------------')
disp(SaveName)
disp('---------------------------------------------------------')