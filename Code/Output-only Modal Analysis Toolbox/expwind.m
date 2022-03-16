function window = expwind(windsize, fs, a)
% window = expwind(windsize, fs) generates an exponential window with a
% given size of - windsize, and a given sampling frequency - fs.
% For your reference:
% The half power point of this window is at t = 0.3464/a.
% The half magnitude point of this wind is at t = 0.6931/a.
%
% syntax: window = expwind(windsize, fs, a)
%
% input: windsize - length of the discrete time window
%        fs - sampling frequency of the recordings
%        a - a real, non-negative number that controls the decay rate  
% output: window - generated exponential window
%
% Hewenxuan Li, March 2021
% Modified on:
% March 31, 2021 - changed window base to the rate of exponential decay
if nargin < 3
    a = 1;
end
if a <= 0
    error('The base to the exponential window function cannot be less than or equal to 0!')
end
dt = 1/fs;
t = linspace(0,windsize*dt,windsize);
window = exp(-a*t);
if size(window,1)<size(window,2)
    window = window';
end

if nargout == 0
    subplot(211)
    plot(window)
    subplot(212)
    [pyy, f] = pwelch(window,boxcar(length(window)),0,[],fs);
    plot(f,10*log10(pyy));
end
end