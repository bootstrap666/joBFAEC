function A = updateTransferMatrix2(sigmaH, Naec, Nbf, A, Nh, M)
% function Anew = updateTransferMatrix2(sigmaH, Naec, Nbf, A, Nh, M)
% random walk model for the LEM plants impulse response
% H = H + deltaH
% where deltaH is a random white gaussian zero mean random matrix (version 2)
% Parameters
%   Anew        -   Updated Transfer Matrix
%   sigmaH      -   variance of deltaH  
%   Naec        -   Acoustic Echo Canceler length
%   Nbf         -   Beamformer Filter length
%   A           -   Current Transfer Matrix
%   Nh          -   Echo channel impulse response length
%   M           -   Number of microphones

A((Naec+M+1):end,2:(Nh+Nbf-1)) = A ((Naec+1):(Naec+M*(Nbf-1)), 1 : (Nh+Nbf-2));

% for i=1:Nbf-1
%     A ((Naec+(i)*M + 1):(Naec+(i+1)*M), (i+1) : (i+Nh)) = A ((Naec+(i-1)*M + 1):(Naec+i*M), i : (i+Nh-1));
% end,

A ((Naec+1):(Naec+M), 1 : Nh) = A ((Naec+1):(Naec+M), 1 : Nh) + sqrt(sigmaH)*randn(M,Nh);
