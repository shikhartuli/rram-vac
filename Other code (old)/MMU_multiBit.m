clc;
close all;
clear;

maxBatchSize = 50;
numBits = 8;
batchSizes = [2:2:maxBatchSize];

mu = 5; % Mean of the switching time distribution
sigma = [1, 2, 4];  % Sigma of the switching time distribution
leftDiscard = 0;
rightDiscard = 4;
offsetMult = 1.1;
numIter = 1000; % 500000;

Clk_M = 1;
Clk_M_Value = 5e-9;
Clk_P = offsetMult * mu * Clk_M; % Clk_P = mu * Clk_M for no accumulation 
% in steady state

dispersion = zeros(length(batchSizes), numIter);

cmap = [1 1 1; hsv];

tic;
for s = 1:1:length(sigma)
for n = 1:numIter
    memoryFill = zeros(maxBatchSize * (2*mu) * Clk_M, numBits);
    memoryFill_base = zeros(maxBatchSize * (2*mu) * Clk_M, numBits);
    memoryTimeDist = zeros(maxBatchSize, numBits);
    last_maxWriteTime_test = 0;
    maxWriteTime_test = 0;
    count = 6;
    for i = 1:1:length(batchSizes)
        count = count - 1;
        last_maxWriteTime_test = last_maxWriteTime_test + maxWriteTime_test;
        % maxWriteTime_test = 0;
        for j = 1:numBits
            writeTime_test = Clk_M * round(mu + sigma(s)*randn);
            while writeTime_test > (mu + rightDiscard * sigma(s)) * Clk_M || writeTime_test < leftDiscard * Clk_M
                writeTime_test = Clk_M * round(mu + sigma(s)*randn);
            end
            if writeTime_test > leftDiscard * Clk_M && writeTime_test < Clk_M
                writeTime_test = Clk_M;
            end
            if i > 1
                memoryTimeDist(i, j) = memoryTimeDist(i-1, j) + writeTime_test;
            else
                memoryTimeDist(i, j) = writeTime_test;
            end
%             if writeTime_test > maxWriteTime_test
%                 maxWriteTime_test = writeTime_test;
%             end
            for k = 1:memoryTimeDist(i, j)
                memoryFill(k, j) = 1 + memoryFill(k, j);
            end
            for k = (last_maxWriteTime_test + 1):1:(last_maxWriteTime_test + writeTime_test)
                memoryFill_base(k, j) = count;
            end
        end
        maxWriteTime_test = 10;
       %         hp = findobj(h, 'type', 'patch');
%         hatch(hp, 45, 'k', '-', 4, 2);
        dispersion(i, n) = sum(abs(memoryTimeDist(i, :) - mean(memoryTimeDist(i, :))))./i;
    end
end
subplot(1, 3, s);
errorbar(batchSizes, mean(dispersion, 2) * Clk_M_Value, (mean(dispersion, 2) - min(dispersion, [], 2)) * Clk_M_Value, (max(dispersion, [], 2) - mean(dispersion, 2)) * Clk_M_Value, '-o');
xlabel('BatchSize');
if s == 1
    ylabel('Normalised dispersion (seconds/word)');
    end
hold on;
plot(batchSizes, mean(dispersion, 2) * Clk_M_Value, 'LineWidth', 2);
fprintf('%0.1f%% Completed\n', s/length(sigma)*100);
hold off;
end
toc;
% subplot(1, 2, 1);
% heatmap(1:numBits, maxBatchSize * (2*mu) * Clk_M:-1:1, flip(memoryFill, 1), 'ColorbarVisible', 'off', 'Colormap', cmap, 'CellLabelColor', 'none');
% subplot(1, 2, 2);
% heatmap(1:numBits, maxBatchSize * (2*mu) * Clk_M:-1:1, flip(memoryFill_base, 1), 'ColorbarVisible', 'off', 'Colormap', cmap, 'CellLabelColor', 'none');


% figure;
% errorbar(batchSizes, mean(dispersion, 2) * Clk_M_Value, (mean(dispersion, 2) - min(dispersion, [], 2)) * Clk_M_Value, (max(dispersion, [], 2) - mean(dispersion, 2)) * Clk_M_Value, '-o');
% hold on;
% errorbar(batchSizes, mean(dispersion, 2), mean(dispersion, 2) - std(dispersion')'./2, mean(dispersion, 2) + std(dispersion')'./2, 'LineStyle', 'None');
% hold off;
% xlabel('BatchSize');
% ylabel('Normalised dispersion (seconds/word)');
