%% Import data from text file
% Script for importing data from the following text file:
%% SELECT DATA SETS
addpath("L:\My Drive\Graduate study\Research\Projects\OS")
[file, path, selectedfile] = openfile('*.mat', 'Select the pre-porcessed modal data file:');
FileName = split(file,'.');
FileName = FileName{1};
load(selectedfile{1})
%% DATA POSTPROCESSING
% Setup Environment
data = impulse_amped_cleansed_data.data;
windsize = impulse_amped_cleansed_data.windsize;
fs = impulse_amped_cleansed_data.fs;

%% Modal FRF esitmates
% ------------------------------------------------------------------------
%                    MODAL IDENFICATION BEGINS
% ------------------------------------------------------------------------
% tfestimate(data_prep(:,8), data_prep(:,2:7), boxcar(windsize),0,[],2000,'Estimator','H1');
InputIndx = 1;
OutputIndx = 2:7;
[FRF, f] = modalfrf(data(:,InputIndx), data(:,OutputIndx), fs, boxcar(windsize),0, 'Sensor', 'dis');

% [FRF, f] = modalfrf(data_prep(:,1), data_prep(:,2:6) - data_prep(:,7), fs, boxcar(windsize),0, 'Sensor', 'dis');

% fn - natural frequencies, dr - damping ratios, ms - mode shapes    
[fn, dr, ms] = modalfit(FRF, f, fs, 1, 'FitMethod', 'lsrf', 'FreqRange', [100 200]);

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

EMA.nf = fn;
EMA.ms = ms;
EMA.dr = dr;
EMA.freq_range = freq_range;
EMA.Nsensor = Nsensors;
EMA.nm = nm;
save([path,'EMA_',FileName,'.mat'],'EMA')
% save([path,FileName,'_cln_amped.mat'],'yamp_test')
%% SAVE DATA TO MAT FORMAT
% ------------------------------------------------------------------------
%            SAVE DATA TO .MAT FORMAT
% ------------------------------------------------------------------------
SaveName = [path, FileName, '.mat'];
save(SaveName, 'data_cln_out')
disp('------------------ Mat file saved ! ---------------------')
disp(SaveName)
disp('---------------------------------------------------------')