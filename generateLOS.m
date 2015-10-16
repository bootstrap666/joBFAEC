function [HLOS, defasagem] = generateLOS(M,N,atrasos, deltaangulos, T0, fs, ...
    PLOS)
% function [HLOS, n0h] = generateLOS(M,N,atrasos, deltaangulos, T0, fs, PLOS)
% Gera a componente devido a linha de visada de um canal de eco sintetico.
% Assume-se inicialmente que o angulo de incidencia da componente devido a
% linha de visada e uma realizacao de uma variavel aleatoria uniforme 
% entre -pi e pi e assume-se que o canal nao e seletivo em frequencia para
% esta componente
% M - Numero de microfones
% N - comprimento do canal
% atrasos - vetor contendo o numero de amostras de defasagem para um angulo
%           de incidencia definido por deltaangulos
% deltaangulos - Indica os intervalos para os quais a defasagem entre um
%                sensor e o proximo e de um numero de amostras
%                correspondente a atrasos(i)
% T0 -  Atraso de transmissao entre a fonte sonora e o microfone mais 
%       proximo a fonte sonora (em segundos)
% fs - Frequencia de amostragem
% PLOS - Potencia do canal de eco na componente de linha de visada
% HLOS - Canal resultante
% n0h - Numero de amostras de atraso de transmissao entre a fonte sonora e
%       o primeiro microfone do array

% sorteia um angulo aleatorio entre 0 e 2*pi uniformemente distribuido
roleta = 2*pi*rand(1);

% sorteia um numero aleatorio -1 ou 1 que define se a componente da linha
% de visada inverte a fase do sinal ou nao
fase = 2*(rand()>0.5)-1;

i=1;
acumulador = 0;
while (acumulador + deltaangulos(i)) < roleta
    acumulador = acumulador + deltaangulos(i);
    i=i+1;
end,

defasagem = atrasos(i);


n0 = ceil(T0*fs);

n0h = max (n0, n0-(M-1)*defasagem);

HLOS = zeros (N,M);

for i=1:M
    HLOS (n0h+(i-1)*defasagem,i) = sqrt(PLOS)*fase;
end,
