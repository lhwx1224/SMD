% ------------------------------------------------------------------------
%                   BEAM SIGNAL AMPLIFICATION EXAMPLE
% ------------------------------------------------------------------------

% Beam Experiments - set up the environment
addpath("L:\My Drive\Graduate study\Research\Projects\OS")
addpath("L:\My Drive\Graduate study\Research\Projects\Output-only Modal Analysis Toolbox\")
% DAQ information
fs = 1000;             % Sampling frequency
fc = 500;              % Cut-off frequency
NCH_cal = 2:7;         % Indices of channels used for calibration from calib data set
NCH_dat = 2:7;         % Indices of channels used for calib and amp from the modal test
RefIndx = 1;           % Reference channel used for calibration (any of the above in NCH_cal)


% Calibrate the data
[data_cal, caldata_cal, rel_mags] = SensorCalibrator(fs, NCH_cal, NCH_dat, RefIndx);

% Assign the calibrated data to a new variable
y = caldata_cal;
y_test = data_cal;
CHref = 8;         % REFERENCE CHANNEL
CHamp = 1:6;       % RESPONSE CHANNELS
% Run amplification
[yamp_test, FRF] = AccAmplifier_TF(y, y_test, CHref, CHamp, fs, fc, 'Nonlinear', 2^10, 1024, [1 400], 'Preamp');

% ------------------------------------------------------------------------
% Compare the soms from the raw and the amplified data
% ------------------------------------------------------------------------

% RAW DATA
[soc, sov, ~, som] = sod(y_test(:,2:7), fs*diff(y_test(:,2:7)));
plot(som./som(1,:))

% AMPED DATA
[soc, sov, ~, som] = sod(yamp_test, fs*diff(yamp_test));
Phi_sod = [som; zeros(1, size(som,1))];
plot(Phi_sod./Phi_sod(1,:))
% ------------------------------------------------------------------------

% Differentiation effect on the harmonic components
[pxx, f] = pwelch(y_test(:,2), boxcar(1024), 512, [], fs);
% [pvv, ~] = pwelch(fs*diff(y_test(:,2)), boxcar(1024), 512, [], fs);
% [paa, ~] = pwelch(fs*diff(fs*diff(y_test(:,2))), boxcar(1024), 512, [], fs);
% [pjj, ~] = pwelch(fs*diff(fs*diff(fs*diff(y_test(:,2)))), boxcar(1024), 512, [], fs);
[pvv, ~] = pwelch(fs*diff(y_test(:,2)), boxcar(1024), 512, [], fs);
[paa, ~] = pwelch(fs*diff(fs*diff(y_test(:,2))), boxcar(1024), 512, [], fs);
[pjj, ~] = pwelch(fs*diff(fs*diff(fs*diff(y_test(:,2)))), boxcar(1024), 512, [], fs);

clf
plot(f, 10*log10(pxx))
hold on
plot(f, 10*log10(pvv))
plot(f, 10*log10(paa))
plot(f, 10*log10(pjj))