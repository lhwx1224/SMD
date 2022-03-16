clear
% Plot stype preoaration
set(0, 'defaultFigureColor', 'w'); 
set(0, 'defaultLineLineWidth', 1.2);
set(0, 'defaultAxesFontName', 'Arial');
set(0, 'defaultTextFontName', 'Arial');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 10);
set(0,'DefaultLineMarkerSize',5);
set(0,'defaultAxesLineWidth',  1);
mycolor = {'#0000b3','#e68a00','#e6b800','#e60000','#009933'};
mycolor1 = {'#003399','#ff9900','#ffcc00','#ff3333','#00e64d'};
% Mycolor 2: 1.Cobalt blue; 2. Love Red; 3. Medium Spring Green; 4. School
% bus yellow; 5. Dimorphotheca Magenta; 6. Dark Violet
mycolor2 = {'#0020C2','#E41B17','#348017','#E8A317','#E3319D','#43C6DB','#842DCE','#0C090A'};
blues = {'#000066','#000099','#0000cc','#0000ff','#3333ff','#6666ff','#9999ff'};
greys = {'#000000','#333333','#666666','#999999','#cccccc','#e6e6e6'};
% Mode Legends
ModeLegends = {'$\widehat{\Phi}_1$','$\widehat{\Phi}_2$','$\widehat{\Phi}_3$', ...
    '$\widehat{\Phi}_4$','$\widehat{\Phi}_5$','$\widehat{\Phi}_6$','$\widehat{\Phi}_7$',...
    '$\widehat{\Phi}_8$','$\widehat{\Phi}_9$','$\widehat{\Phi}_10$'};
%% ------------------ Data Cleansing ----------------------
SavePath = uigetdir('','Set the path to which the figure is saved.');
SavePath = [SavePath,'\'];
[file,~,selectedfile] = openfile('*.mat'); % Obtain file name 
FileName = split(file{:},'.');
FileName = FileName{1};
load(selectedfile{1}); % load data saved from Ace analyzer
fs = 1000;
data = data(:,2:end);
[data_cln_out, windsize] = DataCleansing(data, fs, 1,2:7);
% Save data into a structure which is composed of (1) cleansed data - .data
% (2) propersize of each trial - .windsize and (3) fs - sampling frequency
hammer_data.data = data_cln_out;
hammer_data.windsize = windsize;
hammer_data.fs = fs;
save([SavePath,'cleansed_',FileName,'.mat'],'hammer_data')
%% 
% load('L:\My Drive\Graduate study\Research\My_paper\Journal\Modal_ID_review\Data\cleansed_hammer_test_at_hole_new_config.mat')
% windsize = hammer_data.windsize;
% noverlap = 0;
% freq_range = [0 10; 10 30; 30 100; 100 150];
% data_cln_out = hammer_data.data;
% fs = hammer_data.fs;
% [FRF, f] = modalfrf(data_cln_out(:,1), data_cln_out(:,2:end), fs, boxcar(windsize), noverlap, 'Sensor', 'dis');
% % Peak-picking procedure
% nm = size(freq_range,1);
% fn = zeros(nm,1);
% dr = zeros(nm,1);
% ms = zeros(size(data_cln_out,2)-1,nm);
% for i = 1:nm
%     [fn(i),dr(i),ms(:,i)] = modalfit(FRF,f,fs,1,'FitMethod','lsce','FreqRange',freq_range(i,:));
% end
%% Modal Isolation and Identification
% disp('Select the Data Set Need to be Adjusted:')
[file, path, DataName] = openfile('*.mat','- select the modal testing data:');
load(DataName{1})
windsize = impulse_cleansed_data.windsize;
fs = impulse_cleansed_data.fs;
data_cln_out = impulse_cleansed_data.data;
noverlap = 0;
% 
% freq_range = [0 100; 50 200; 100 450; 300 500];
freq_range = [0 40; 30 120];
[FRF, f] = modalfrf(data_cln_out(:,1), data_cln_out(:,2:end), fs, boxcar(windsize), noverlap, 'Sensor', 'dis');
% Peak-picking procedure
nm = size(freq_range,1);
fn = zeros(nm,1);
dr = zeros(nm,1);
ms = zeros(size(data_cln_out,2)-1,nm);
for i = 1:nm
    [fn(i),dr(i),ms(:,i)] = modalfit(FRF,f,fs,1,'FitMethod','lsce','FreqRange',freq_range(i,:));
end
% modalfit(FRF,f,fs,6);
%% Multi-mode identification
[file, path, DataName] = openfile('*.mat','- select the modal testing data:');
load(DataName{1})
windsize = impulse_cleansed_data.windsize;
fs = impulse_cleansed_data.fs;
data_cln_out = impulse_cleansed_data.data;
noverlap = 0;

[FRF, f] = modalfrf(data_cln_out(:,1), data_cln_out(:,2:end), fs, boxcar(windsize), noverlap, 'Sensor', 'dis');
[fn, dr, ms] = modalfit(FRF, f, fs, 7, 'FitMethod', 'lsrf', 'FreqRange', [0 800]);
% [fn, dr, ms] = modalfit(FRF, f, fs, 6, 'FitMethod', 'lsrf', 'FreqRange', [0 800]);

figure(2),clf
modalsd(FRF,f,fs,'MaxModes',30);

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
legend(p1,ModeLegends,'Location','EastOutside','NumColumns',2)
title('Cleansed-Only Data Set')
%% Multi-mode identification
addpath("L:\My Drive\Graduate study\Research\Projects\OS")
[file, path, DataName] = openfile('*.mat','- select the modal testing data:');
load(DataName{1})
windsize = impulse_amped_cleansed_data.windsize;
fs = impulse_amped_cleansed_data.fs;
data_cln_out = impulse_amped_cleansed_data.data;
noverlap = 0;

DataInspection(data_cln_out)

[FRF, f] = modalfrf(data_cln_out(:,1), data_cln_out(:,2:end), fs, boxcar(windsize), noverlap, 'Sensor', 'dis');
[fn, dr, ms] = modalfit(FRF, f, fs, 4, 'FitMethod', 'lsrf', 'FreqRange', [0 200]);
[fn, dr, ms] = modalfit(FRF, f, fs, 2, 'FitMethod', 'lsrf', 'FreqRange', [400 800]);
[fn, dr, ms] = modalfit(FRF, f, fs, 7, 'FitMethod', 'lsrf', 'FreqRange', [0 800]);

figure(2),clf
modalsd(FRF,f,fs,'MaxModes',30);

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
legend(p1,ModeLegends,'Location','EastOutside','NumColumns',2)
title('Amplified-and-Cleansed Data Set')
%% Visualization of the mode shape and the modal parameters
Phi = real(ms);
Phi_norm = Phi./Phi(1,:);
Phi_norm = [Phi_norm; zeros(1,size(Phi_norm,2))];
figure(1),clf
p1 = plot(Phi_norm);
hold on
plot([1 size(Phi_norm,1)],[0 0],'--','Color','#B6B6B2')
pbaspect([2 1 1])
title('Experimental Modal Analysis --- LSCE')
xlabel('Sensor Location')
ylabel('Normalized Magnitude')
grid on
legend(p1,ModeLegends,'Location','SouthWest')

% SavePath = uigetdir('','Set the path to which the figure is saved.');
% SavePath = [SavePath,'\'];
% SaveName = 'EMA_LSCE_mar_9_impulse3_4db';
% exportgraphics(figure(1),[SavePath,[SaveName,'.png']],'Resolution',600)

plot(data_cln_out)
title('Cleansed Data')
xlabel('Sample')
ylabel('Magnitude (V)')
pbaspect([2 1 1])
axis tight
% SaveName = 'CleansedData_mar_9_impulse3_4db';
% exportgraphics(figure(1),[SavePath,[SaveName,'.png']],'Resolution',600)

[pyy, x] = pwelch(data_cln_out,boxcar(windsize), noverlap,[],fs);
plot(x, 10*log10(pyy))
pbaspect([2 1 1])
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
title('Power Spectral Density Estimate')
% SaveName = 'PSD_mar_9_impulse3_4db';
% exportgraphics(figure(1),[SavePath,[SaveName,'.png']],'Resolution',600)

H = pyy(:,1:end-1)./pyy(:,end);
plot(x,H)
set(gca,'yscale','log')
pbaspect([2 1 1])
xlabel('Frequency (Hz)')
ylabel('FRF (dB Magnitude)')
title('Frequency Response Function Estimate')
grid on
% SaveName = 'FRF_mar_9_impulse3_4db';
% exportgraphics(figure(1),[SavePath,[SaveName,'.png']],'Resolution',600)

% modalfrf(data_cln(:,end), data_cln(:,1:end-1), fs, boxcar(windsize), noverlap, 'Sensor', 'dis');
% figure(1)
% pbaspect([4 1 1])
% SaveName = 'FRF_mar_9_impulse3_4db';
% exportgraphics(figure(1),[SavePath,[SaveName,'.png']],'Resolution',600)