function A = transferMatrix(H, Naec, Nbf)
% function A = transferMatrix(H, Naec, Nbf)
% Generates the transfer matrix of a multiple microphone acoustic echo
% channel 
% H - Echo channel matrix
% Naec - echo canceller length
% Nbf - beamformer length

M = size(H,2);
Nh = size(H,1);

A = zeros (Naec + M*Nbf, Nh+Nbf + M*Nbf-1);

A(1:Naec, 1:Naec) = -eye(Naec);

for i=1:Nbf
    A ((Naec+(i-1)*M + 1):(Naec+i*M), i : (i+Nh-1)) = H';
end,

A((Naec+1):(Naec + M*Nbf),(Nh + Nbf):(Nh + Nbf + M*Nbf-1)) ...
    = eye(M*Nbf);