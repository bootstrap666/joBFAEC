function fig = plotsimresults(DMonteCarlo, Dstatmodel, Jinf, iterations, decimationFactor)

if (nargin < 5)
    decimationFactor = 1;
end
if (nargin < 4)
    iterations = length(DMonteCarlo);
end


fig = figure();
set(gca,'FontSize',14);
plot(1:decimationFactor:iterations, 10*log10(DMonteCarlo(1:decimationFactor:iterations)));
hold on;
plot((1:length(Dstatmodel))*iterations/length(Dstatmodel), 10*log10(Dstatmodel),'r','LineWidth',2);
plot([1 iterations], 10*log10(Jinf)*[1 1],'r:','LineWidth',2);
xlim([1 iterations]);
xlabel('n', 'FontSize',14);
ylabel('E\{d^2[n]\}', 'FontSize',14);
grid on;
legend('MC simulation','model prediction', 'steady-state MOP');