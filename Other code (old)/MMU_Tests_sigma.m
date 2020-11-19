clc;
close all;
clear;

P_cycles = 20000;

mu = 20;
sigma_min = 1;
sigma_max = 10;
sigmaArray = (sigma_min:1:sigma_max);
offsetMult = 1.1;
numIter = 1;
maxBufferArray = zeros(length(sigmaArray), numIter);
BufferWithTimeArray = {};
writeTimes = {};
leftDiscard = 0;
rightDiscard = 3;


disp('Simulation Running');
for i = 1:1:numIter
    for j = 1:1:length(sigmaArray)
        sigma = sigmaArray(j);
        output = MMU(P_cycles, mu, sigma, leftDiscard, rightDiscard, offsetMult);
        maxBufferArray(j, i) = output{1};
        BufferWithTimeArray{j} = output{2}(:, [1, 2]);
        writeTimes{j} = output{3};
    end
    % fprintf('%d%% Completed\n', i/numIter*100);
end

fig1 = figure;
errorbar(sigmaArray, mean(maxBufferArray, 2), mean(maxBufferArray, 2) - min(maxBufferArray, [], 2), max(maxBufferArray, [], 2) - mean(maxBufferArray, 2), '-o')
xlabel('\sigma')
ylabel('Required Buffer Size')
saveas(fig1, sprintf('Req_BufferSize_with_sigma_mu%d_offsetMult%0.1f.fig', mu, offsetMult));
saveas(fig1, sprintf('Req_BufferSize_with_sigma_mu%d_offsetMult%0.1f.png', mu, offsetMult));

fig2 = figure;
for i = 1:1:length(sigmaArray)
    plot(BufferWithTimeArray{i}(:, 1), BufferWithTimeArray{i}(:, 2), 'DisplayName', sprintf('sigma = %d', sigmaArray(i)), 'LineStyle', '-.')
    hold on
end
xlabel('Time (Clk_{M}s)')
ylabel('Buffer Size')
legend show
grid on;
saveas(fig2, sprintf('BufferSize_with_time_(different_sigma)_mu%d_offsetMult%0.1f.fig', mu, offsetMult));
saveas(fig2, sprintf('BufferSize_with_time_(different_sigma)_mu%d_offsetMult%0.1f.png', mu, offsetMult));

fig3 = figure;
for i = 1:1:length(sigmaArray)
    histogram(writeTimes{i}, 'DisplayName', sprintf('sigma = %d', sigmaArray(i)), 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3)
    hold on
end
xlabel('Time (Clk_{M}s)')
ylabel('Count')
legend show
saveas(fig3, sprintf('Histograms_(different_sigma)_mu%d_offsetMult%0.1f.fig', mu, offsetMult));
saveas(fig3, sprintf('Histograms_(different_sigma)_mu%d_offsetMult%0.1f.png', mu, offsetMult));
