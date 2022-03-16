%% Modal testing data pre-processing program for sine-swept (SS)
% This program is used for pre-processing collected modal testing data from
% sine-swept excitations. The following steps will be followed:
% (1) Data selection
% (2) Data calibration (mean-std adjustment according to collocation tests)
% (3) Data amplification (adjust the energy bias/nonlinearity in the
%     collected data according to coloocation tests from calibration data 
%     under wide-band excitation)   
% (4) Fast preanalysis 
% (5) Data saving
%% PART I: SELECT DATA SETS (THIS PART IS SUBJECTED TO CHANGE)
% addpath("L:\My Drive\Graduate study\Research\Projects\OS")
% [file, path, selectedfile] = openfile('*.txt', 'Select the modal data file:');
% FileName = split(file,'.');
% FileName = FileName{1};
% % Set up the Import Options and import the data
% opts = delimitedTextImportOptions("NumVariables", 9);
% 
% % Specify range and delimiter
% opts.DataLines = [1, Inf];
% opts.Delimiter = "\t";
% 
% % Specify column names and types
% opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"];
% opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double"];
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Import and concatenate the data
% data = cell(length(selectedfile),1);
% for i = 1:length(selectedfile)
%     data{i} = table2array(readtable(selectedfile{i}, opts));
% end
% 
% data = cell2mat(data);
%% DATA PREPROCESSING
% ------------------------------------------------------------------------
%                       PART II - DATA CALIBRATION
% ------------------------------------------------------------------------
% Setup Environment
addpath("L:\My Drive\Graduate study\Research\Projects\OS")
Amp = 0;        % Use amplification?
fs = 2000;      % Sampling frequency
fc = 1007;      % Cutoff frequency
windsize = 1024; % Window size for spectral density estimation
NCH_cal = 1:7;  % Channel number for the calibration data set
NCH_dat = 1:7;  % Channel number for the modal test data set
RefIndx = 1;    % Reference channel used in the calibration data set

[data_cal, caldata_cal, rel_mags] = SensorCalibrator(fs, NCH_cal, NCH_dat, RefIndx);

% ------------------------------------------------------------------------
%                     PART III - Run amplification
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
% prepared data as the amplified data
data_prep = yamp_test;
else
    data_prep = data_cal; % IF not amplified, use the calibrated data directly
end

%
spectrogram(data_prep(:,8), boxcar(1024), 512, [], 2000)
mycolorbar()
pbaspect([1.618 1 1])
set(gcf, 'Papersize', [6 6*0.5])
set(gcf, 'PaperPosition', [0 0 6 6*0.5])
% 
tfestimate(data_prep(:,1), data_prep(:,2:7), boxcar(windsize),0,[],2000,'Estimator','H1');
set(gcf, 'Papersize', [6 6*0.5])
set(gcf, 'PaperPosition', [0 0 6 6*0.5])

%% Modal FRF esitmates
% ------------------------------------------------------------------------
%                 PART IV - MODAL IDENFICATION BEGINS
% ------------------------------------------------------------------------
% tfestimate(data_prep(:,8), data_prep(:,[2:7, 9]), boxcar(windsize),0,[],2000,'Estimator','H1');
% tfestimate(data_prep(:,9), data_prep(:,2:7), boxcar(windsize),0,[],2000,'Estimator','H1');

[FRF, f] = modalfrf(data_prep(:,9), data_prep(:,[2:7]), fs, boxcar(windsize),0, 'Sensor', 'dis');
% [FRF, f] = modalfrf(data_prep(:,8), data_prep(:,[2:7]) - data_prep(:,9), fs, boxcar(windsize),0, 'Sensor', 'dis');
% [FRF, f] = modalfrf(data_prep(:,1), data_prep(:,2:6) - data_prep(:,7), fs, boxcar(windsize),0, 'Sensor', 'dis');

% fn - natural frequencies, dr - damping ratios, ms - mode shapes
[fn, dr, ms] = modalfit(FRF, f, fs, 1, 'FitMethod', 'lsrf', 'FreqRange', [150 250]);

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

% save([path,FileName,'_cln_amped.mat'],'yamp_test')
%% SAVE DATA TO MAT FORMAT
% ------------------------------------------------------------------------
%            PART V - SAVE DATA TO .MAT FORMAT
% ------------------------------------------------------------------------
SaveName = [path, FileName, '.mat'];
save(SaveName, 'data_cln_out')
disp('------------------ Mat file saved ! ---------------------')
disp(SaveName)
disp('---------------------------------------------------------')