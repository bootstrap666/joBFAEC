function [C, f] = generateRestrMatrixAng (angle, desiredimpres, dmic, c,...
    fs, M, threshold)
% function [C, f] = generateRestrMatrixAng (angle, desiredimpres, dmic, c,...
%    fs, M, limiar)
% Generates the constraint matrix and constraint frequency response of an 
% uniform linear microphone array using the method described in Buckley, 
% 1987, Spatial/Spectral filtering with linearly constrained minimum 
% variance beamformers, IEEE Transactions on Acoustics, Speech and  Signal 
% Processing, v 35, n 3, pp 249 - 266.
%
% angle - Desired angle of arrival (in radians)
% desiredimpres - Desired impulse response in the angle of arrival 
% dmic - distance between adjacent microphones (in meters)
% c - speed of sound (typically 343 m/s)
% fs - sampling rate (in samples/sec)
% M - number of microphones
% threshold - minimal output power approximation in the least squares
%       process (default 0.99)
% C - constraint matrix
% f - transformed constraint frequency response where
% C'*w = f 

if (nargin<7)
    threshold=0.99;
end,
    
N = length(desiredimpres);

[respfreqdes, omega] = freqz(desiredimpres,1,N) ;

tau = (((0:M-1) - (M-1)/2)*dmic*sin(angle)*fs/c)';

theta = zeros (M*N,N);

for i=1:N
    theta((i-1)*M+1:i*M,:) = (tau-i+1)*omega';
end,

% Gera a matriz de restricoes com 2*N restricoes reais
Atheta = [cos(theta), sin(theta)];

% Separa a parte real e imaginaria da resposta em frequencia desejada para
% obter restricoes reais

F = [real(respfreqdes); imag(respfreqdes)];

[U, S, V] = svd(Atheta,0);

% Encontra o numero de valores singulares necessarios para obter
% limiar*pottotal de potencia para definir o tamanho do subespaco de sinal
pottotal = trace(S);

nSV=0;
sum = 0.0;
while (sum<threshold*pottotal)
    nSV = nSV+1;
    sum = sum + S(nSV,nSV);
    if (nSV == N)
        break;
    end,
end,

% Separa somente a projecao das restricoes no subespaco de sinal

Sl = S(1:nSV,1:nSV);
Vl = V(:,1:nSV);
Ul = U(:,1:nSV);

C = Ul;
f = (Sl\Vl')*F;

