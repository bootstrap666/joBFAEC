function H = generateindepchannels (M, n0, TR60, fs, power, N)
% function H = generateindepchannels (M, T0, TR60, fs, power, N)
% Generates M synthetic random independent echo channels from a sequence of
% independent white Gaussian random variables. Assums the exponential decay 
% model of eq. 7 in Breining et al 1999
% Acoustic Echo Control An application of Very-High-Order Adaptive Filters
% - IEEE Signal processing magazine pp 42--69
% Parameters
%   M       -   Number of microphones
%   n0      -   Samples delay from the loudspeaker to the center of the microphone
%               array. It can be computed as n0 = fs*d/c where d is the
%               loudspeaker-microphone distance, c is the speed of sound
%               (about 343 m/s).
%   TR60    -   Reverberation time
%   fs      -   Sampling rate (default=8192)
%   power   -   power gain of each echo channel (default=ones(M,1))
%   N       -   Echo impulse response length(default=2^nextpow2(n0 + ceil(TR60*fs)))
%   H       -   Resulting echo channel

if nargin < 3
    fs = 2^13;
end,

if nargin < 4
    power = ones(M,1);
end,

if nargin < 5
    N = 2^nextpow2(n0 + ceil(TR60*fs));
end,

% Time constant
alpha = -3/(fs*TR60*log10(exp(1)));

Neff = N-(n0-1);

sigmah = power/sum(exp(2*alpha*(0:(Neff-1))'));

H = zeros (N,M);
R = randn(Neff,M);

for i=1:M
    H(:,i) = [zeros(n0-1,1); exp(alpha*(1:Neff)') .* (sqrt(sigmah(i))*R(:,i))];
end,
