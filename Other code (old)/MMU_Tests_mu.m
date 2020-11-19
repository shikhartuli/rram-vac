clc;
close all;
clear;

P_cycles = 20000;

mu = 
muArray = mu_min:2:mu_max;
sigma = 5;
offsetMult = 1.1;
numIter = 1;
maxBufferArray = zeros(length(muArray), numIter);
BufferWithTimeArray = {};
writeTimes = {};
leftDiscard = 0;
rightDiscard = 3;


disp('Simulation Running');
for i = 1:1:numIter
    parfor j = 1:1:length(muArray)
        mu = muArray(j);
        output = MMU(P_cycles, mu, sigma, leftDiscard, rightDiscard, offsetMult);
        maxBufferArray(j, i) = output{1};
        BufferWithTimeArray{j} = output{2}(:, [1, 2]);
        writeTimes{j} = output{3};
        disp(j);
    end
    fprintf('%d%% Completed\n', i/numIter*100);
end

fig1 = figure;
errorbar(muArray, mean(maxBufferArray, 2), mean(maxBufferArray, 2) - min(maxBufferArray, [], 2), max(maxBufferArray, [], 2) - mean(maxBufferArray, 2), '-o')
xlabel('\mu')
ylabel('Required Buffer Size')
saveas(fig1, sprintf('Req_BufferSize_with_mu_sigma%d_offsetMult%0.1f.fig', sigma, offsetMult));
saveas(fig1, sprintf('Req_BufferSize_with_mu_sigma%d_offsetMult%0.1f.png', sigma, offsetMult));

fig2 = figure;
for i = 1:1:length(muArray)
    plot(BufferWithTimeArray{i}(:, 1), BufferWithTimeArray{i}(:, 2), 'DisplayName', sprintf('mu = %d', muArray(i)))
    hold on
end
xlabel('Time (Clk_{M}s)')
ylabel('Buffer Size')
legend('NumColumns', 2)
legend show
grid on;
saveas(fig2, sprintf('BufferSize_with_time_(different_mu)_sigma%d_offsetMult%0.1f.fig', sigma, offsetMult));
saveas(fig2, sprintf('BufferSize_with_time_(different_mu)_sigma%d_offsetMult%0.1f.png', sigma, offsetMult));

fig3 = figure;
for i = 1:1:length(muArray)
    histogram(writeTimes{i}, 'DisplayName', sprintf('mu = %d', muArray(i)), 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3)
    hold on
end
xlabel('Time (Clk_{M}s)')
ylabel('Count')
legend('NumColumns', 2)
legend show
saveas(fig3, sprintf('Histograms_(different_mu)_sigma%d_offsetMult%0.1f.fig', sigma, offsetMult));
saveas(fig3, sprintf('Histograms_(different_mu)_sigma%d_offsetMult%0.1f.png', sigma, offsetMult));