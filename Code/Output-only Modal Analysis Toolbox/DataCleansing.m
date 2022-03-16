function [data_cln_out, windsize] = DataCleansing(data, fs, input_indx, output_indx, A, windratio, FullOutputData, VarianceThreshold)
% [data_cln, windsize] = DataCleansing(data, fs, input_indx, output_indx,
% A, windratio) returns the cleansed data from raw field measurement data
% whose columns represent measurement in time. Given the variable data, the
% input is designated using the parameter input_indx, the output channels
% are chosen by output_indx; as for windowing, an exponential window decay
% rate is determined by A and the windowed data size is determined by the
% variable, windratio.
%
% syntax: [data_cln, windsize] = DataCleansing(data, fs, input_indx,
% output_indx, A, windratio) 
%
% input: data - raw data with each column stores measurement of each
%        channel
%        fs - sampling frequency used in the experiment
%        input_indx - 1d array stores indices of the input channels
%        output_indx - 1d array stores indices of the output channels,
%        e.g., output_indx = [2:5] -> if column 2 to 5 are output data
%        A - rate of exponential decay used in the exponential window
%        windratio - a real number between 0 and 1 indicating the ratio of
%        the truncated-and-cleansed data remained and forced to zero at its
%        tail.
% output: data_cln_out - cleansed data set
%         windsize - length of the cleansed-and-truncated data in between
%         two impulses.
%
% Copyright by Hewenxuan Li Mar 31, 2021
% As a part of the output-only modal analysis toolbox
% Last Modified:
% Oct. 4, 2021 - Assigned full cleansed data set as the output variable
% data_cln_out on Line 155.
% Oct. 7, 2021 - Added FullOutputData output option, if a fullsize data is
% needed, one can input FullOutputData = 'Full' as input to the function.

% SETUP THE ENVIRONMENT 
set(0,'DefaultFigureWindowStyle','docked')
addpath("L:\My Drive\Graduate study\Research\Projects\OS")
mystyle();

% ----- Modified 10/7/2021 ------
if nargin < 8
    FullOutputData = 'False';
    VarianceThreshold = 100;
elseif nargin < 6
    windratio = 0.5;
    FullOutputData = 'False';
    VarianceThreshold = 100;
elseif nargin < 5
    A = 1;
    windratio = 0.5;
    FullOutputData = 'False';
    VarianceThreshold = 100;
end

% -------------------------------

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

% ------------------ Data Cleansing ----------------------
% If condition on whether a time vector exists
if round(1/mean(data(2:11,1) - data(1:10,1))) == fs
    time = data(:,1);
    data(:,1) = [];
    disp('Time vector in calibration data detected and removed for amplification!')
    time_flag = 1;
    impulse = data(:, input_indx);
    data_out = data(:,output_indx);
else
    time_flag = 0;
    disp('No time vectors identified.')
    impulse = data(:, input_indx);
    data_out = data(:,output_indx);
end
[~, n] = size(data);

% High-pass filtering with Butterworth filter
fc = 3;
[b,a] = butter(6, fc/(fs/2), 'high');
impulse_hp = filtfilt(b,a,impulse);
all_indx = 1:1:n;

% excluded channels for amplification (assuming the first channel is the time vector)
if time_flag
    CH_exclude = ~ismember(all_indx, [1 input_indx output_indx]); 
else
    CH_exclude = ~ismember(all_indx, [input_indx output_indx]); 
end
CH_include = ismember(all_indx, [input_indx output_indx]); % included channels for amplification
data_exclude = data(:, CH_exclude);

% -------------------------------------------------------------------------
% Visualization for high-pass filtering the input signal (debugging mode) -
% plot(impulse)
% hold on
% plot(impulse_hp)
% -------------------------------------------------------------------------

% Impulse detection
indx_imp = impulse_hp>VarianceThreshold*std(impulse_hp(1:1000)); % Variance based detection
plot(impulse_hp)
hold on
plot(indx_imp)
plot(data(:,2))
% Extract the impact location based on the relative variance between
% impulses and the noise.
loc_imp = [];
for i = 1:length(indx_imp)
    if indx_imp(i) == 1 && indx_imp(i-1) == 0
        loc_imp = [loc_imp i];
    end
end
loc_imp = [loc_imp size(data,1)]; % Location of the impulse
loc_imp = loc_imp - 1;
% Sift out the improper impulses
dloc = loc_imp(2:end) - loc_imp(1:end-1); % Time interval between imp.s
flag = [];
for i = 1:length(dloc)
    if dloc(i) < 2*fs
        flag = [flag, i + 1];
    end
end
loc_imp(flag) = [];
% Extract the data sets from the location indices
y = cell(length(loc_imp)-1,1);
min_size = [];
for i = 1:length(loc_imp)-1
    % Changed impulse to impulse_hp 10/7/21 for the following line
    y{i} = [impulse(loc_imp(i):loc_imp(i+1)), data_out(loc_imp(i):loc_imp(i+1),1:end)];
    if i == 1
        min_size = length(y{i});
    else 
        if length(y{i})<min_size
            min_size = length(y{i});
        end
    end
end

% Application of the exponential window
windsize = fix(min_size*windratio);
wind = expwind(windsize,fs,A);
% Visualize the selected window
figure(1),clf
subplot(311)
plot(wind)
grid on
axis tight
xlabel('Sample')
ylabel('Magnitude')
title('Window for Data Cleansing')
subplot(312)
pwelch(wind,boxcar(size(wind,1)),0)
% Matrix window -> fit the data size
wind = repmat(wind,1,size(y{1},2));
wind(:,1) = ones(size(wind,1),1);
% Data inspection for all impulse events
subplot(313)
plot(y{1}(1:windsize,end),'color','k');
hold on
plot(wind.*y{1}(1:windsize,end),'color','r');
axis tight
xlabel('Sample')
ylabel('Magnitude')
title('Windowing Effect')
legend('Raw-and-truncated data','Windowed-and-truncated data')
% Truncate the data and put them back into a matrix
data_cln = zeros(windsize*length(y),size(y{1},2));
for i = 1:length(y)
    data_cln((i-1)*windsize+1:i*windsize,1:end) = wind.*y{i}(1:windsize,1:end);
%     data_cln((i-1)*length(wind)+1:i*length(wind),end) = y{i}(1:length(wind),end);
end
figure(3),clf
subplot(311)
plot(data_cln(:,2:end))
ylabel('Magnitude')
xlabel('Sample')
title('Cleansed Response Data')
axis tight
subplot(312)
plot(data_cln(:,1))
ylabel('Magnitude')
xlabel('Sample')
title('Cleansed Input Data')
axis tight
subplot(313)
pwelch(data_cln,boxcar(windsize),0,[],fs);
retain_indx = 1:1:size(data_cln,2);

while_count = 0;
while ~isequal(lower(retain_indx),'all')
    disp(['========== There are ', num2str(size(loc_imp,2)), ' impulses identified! =========='])
    retain_indx = input('Please select the desirable impulses in the cleansed \n time series (e.g., 1:5, if all entries are desirable, enter string "All"):');
    if isequal(lower(retain_indx), "all")
        if while_count == 0
        disp('Output impulse responses done! All impulses will be saved!')
        data_cln_out = data_cln;
        break
        else 
            disp(['Output the following impulse response histories:', num2str(retain_indx)])
        break
        end
    else
        data_cln_out = [];
    end
    for i = 1:length(retain_indx)
        data_cln_out = [data_cln_out; data_cln((retain_indx(i)-1)*windsize+1:(retain_indx(i))*windsize,:)];
    end
    figure(3),clf
    subplot(211)
    plot(data_cln_out)
    ylabel('Magnitude')
    xlabel('Sample')
    title('Cleansed Data Ready for Output')
    axis tight
    subplot(212)
    pwelch(data_cln_out,boxcar(windsize),0,[],fs);
    while_count = while_count + 1;
    retain_indx_disp = retain_indx; % Save for display purpose
end
sgtitle('Data Cleansing Results')

% ----- Modified 10/7/2021 ------ WORK IN PROGRESS!
if isequal(lower(FullOutputData), 'full')
    warning('NOTE: the input is located on the second column in the cleansed data!')
    disp('Full data will be returned!')
    if time_flag
        newtime = time(1:size(data_cln_out,1));
        data_out_full = zeros(size(data_cln_out,1),n);
    else
        newtime = linspace(0, 1/fs*(size(data_cln_out,1)-1), size(data_cln_out,1));
        data_out_full = zeros(size(data_cln_out,1),n + 1);
    end
    data_out_full(:,1) = newtime;
    data_out_full(:,CH_include) = data_cln_out;
    data_out_full(:,CH_exclude) = data_exclude(1:size(data_cln_out,1),:);
    data_cln_out = data_out_full;
else
    warning('NOTE: the input is located on the first column in the cleansed data!')
end
% -------------------------------
