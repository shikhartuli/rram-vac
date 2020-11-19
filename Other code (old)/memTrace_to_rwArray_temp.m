clc;
% close all;
% clear all;


mu = 5;
sigma = 1; 
leftDiscard = 0;
rightDiscard = 4;
numBits = 32;
Clk_M = 1;
Clk_M_Value = 5e-9; % in seconds
bitReadTime = 1* Clk_M;
switchDetectTime = 0.2 * Clk_M; % 1ns
writeVolt = 1; % in Volts
I_LRS = 100e-6; % in Amps
I_HRS = 15e-6; % in Amps
readEnergyPerBit = 1e-12; % in J
bufferAccessEnergy = 400e-15;
debug = 0;

writeBufferSize = 10;
waitBufferSize = 80;

numIter =  5;
offsetMult = 1; % 1:0.1:2;
offsetMultBase = 2; % (mu + rightDiscard * sigma) / mu;
writeStallTimes_parallel = zeros(numIter, 1);
readStallTimes_parallel = zeros(numIter, 1);
times_parallel = zeros(numIter, 1);
readEnergy_parallel = zeros(numIter, 1);
writeEnergy_parallel = zeros(numIter, 1);
readStallTimes_parallel_base = zeros(numIter, 1);
writeStallTimes_parallel_base = zeros(numIter, 1);
times_parallel_base = zeros(numIter, 1);
readEnergy_parallel_base = zeros(numIter, 1);
writeEnergy_parallel_base = zeros(numIter, 1);

% timesStack = [];
% energyStack = [];
timesStack_norm = [];
energyStack_norm = [];

groupLabels = { 'CS', 'CS(32)', 'FE', 'DT', 'DT(C)', 'FE+DT', 'MM', 'Conv' };

% memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportCS_full.csv');
% rwArray = {};
% % Forming rwArray from memTrace
% cellArray = textscan(memTrace, '%c %d %s');
% rwArray = cell(length(cellArray{1, 1}), 2);
% for k = 1:1:length(cellArray{1, 1})
%     rwArray{k, 1} = cellArray{1, 1}(k);
%     if cellArray{1, 1}(k) == 'w'
%         rwArray{k, 2} = [cellArray{1, 2}(k), hex2dec(cellArray{1, 3}(k))];
%     else
%         rwArray{k, 2} = cellArray{1, 2}(k);
%     end
% end

tic;
for exp = 6
    if exp == 1
        memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportCS_full.csv');
        disp('Runnning CS');
    elseif exp == 2
        memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportCS_full_2.csv');
        disp('Runnning CS (2)');
    elseif exp == 3
        memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportCS_full_4.csv');
        disp('Runnning CS (4)');
    elseif exp == 4
        memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportCS_full_8.csv');
        disp('Runnning CS (8)');
    elseif exp == 5
        memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportCS_full_16.csv');
        disp('Runnning CS (16)');
    elseif exp == 6
        memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportCS_full_32.csv');
        disp('Runnning CS (32)');
        
%     elseif exp == 2
%         memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportFE.csv');
%         disp('Running FE');
%     elseif exp == 3
%         memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportDT.csv');
%         disp('Running DT');
%     elseif exp == 4
%         memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportDT_cached_2.csv');
%         disp('Running DT (Cached)');
%     elseif exp == 5
%         memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportFE+DT.csv');
%         disp('Running FE+DT');
%     elseif exp == 6
%         memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportMM_new.csv');
%         disp('Running MM');
%     elseif exp == 7
%         memTrace = fopen('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 10 Apps results\\reportConv.csv');
%         disp('Running Conv');
    
    end
    rwArray = {};
    % Forming rwArray from memTrace
    cellArray = textscan(memTrace, '%c %d %s');
    rwArray = cell(length(cellArray{1, 1}), 2);
    for k = 1:1:length(cellArray{1, 1})
        rwArray{k, 1} = cellArray{1, 1}(k);
        if cellArray{1, 1}(k) == 'w'
            rwArray{k, 2} = [cellArray{1, 2}(k), hex2dec(cellArray{1, 3}(k))];
        else
            rwArray{k, 2} = cellArray{1, 2}(k);
        end
    end
    parfor i = 1:numIter
    %     memTrace = fopen(sprintf('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 6 MMU CS\\ECG_3s_each\\report_global_%d.csv', randi([0, 19])));
    %     memTrace = fopen(sprintf('/scrap/users/tuli/ECG_3s_each/report_global_%d.csv', randi([0, 19])));
    %     rwArray = {};
    %     % Forming rwArray from memTrace
    %     cellArray = textscan(memTrace, '%c %d %s');
    %     rwArray = cell(1000, 2); % cell(length(cellArray{1, 1}), 2);
    %     for k = 1:1000 % length(cellArray{1, 1})
    %         rwArray{k, 1} = cellArray{1, 1}(k);
    %         if cellArray{1, 1}(k) == 'w'
    %             rwArray{k, 2} = [cellArray{1, 2}(k), hex2dec(cellArray{1, 3}(k))];
    %         else
    %             rwArray{k, 2} = cellArray{1, 2}(k);
    %         end
    %     end
        [writeStallTimes_parallel(i), readStallTimes_parallel(i), times_parallel(i), readEnergy_parallel(i), writeEnergy_parallel(i)] = MMU_rw_parallel_wEnergy(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult, writeBufferSize, waitBufferSize, bitReadTime, switchDetectTime, writeVolt, I_LRS, I_HRS, readEnergyPerBit, bufferAccessEnergy, Clk_M_Value, 0, debug);
        [writeStallTimes_parallel_base(i), readStallTimes_parallel_base(i), times_parallel_base(i), readEnergy_parallel_base(i), writeEnergy_parallel_base(i)] = MMU_rw_parallel_wEnergy(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMultBase, 1, waitBufferSize, bitReadTime, switchDetectTime, writeVolt, I_LRS, I_HRS, readEnergyPerBit, bufferAccessEnergy, Clk_M_Value, 1, debug);
    %     [times_serial(1, j), readEnergy_serial(i, j), writeEnergy_serial(i, j), maxWaitBufferSize_serial(i, j)] = MMU_rw_serial_wEnergy(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult * numBits, writeBufferSize(i), bitReadTime, switchDetectTime, writeVolt, I_LRS, I_HRS, readEnergyPerBit, bufferAccessEnergy, Clk_M_Value, 0, debug);
    end
    timesStack(exp, 1, 1) = mean(times_parallel, 1) * Clk_M_Value;
    timesStack(exp, 2, 1) = mean(times_parallel_base, 1) * Clk_M_Value;
    energyStack(exp, 1, 1) = mean(readEnergy_parallel, 1);
    energyStack(exp, 2, 1) = mean(readEnergy_parallel_base, 1);
    energyStack(exp, 1, 2) = mean(writeEnergy_parallel, 1);
    energyStack(exp, 2, 2) = mean(writeEnergy_parallel_base, 1);
    toc;
end

for exp = 1:1:8
    timesStack_norm(exp, :) = timesStack(exp, :)/(timesStack(exp, 2));
    energyStack_norm(exp, :, :) = energyStack(exp, :, :)/(energyStack(exp, 2, 1) + energyStack(exp, 2, 2));
end

plotBarStackGroups(timesStack, groupLabels);
xlabel('Applications');
ylabel('Total Process Time (seconds)');
box on;

plotBarStackGroups(energyStack, groupLabels);
xlabel('Applications');
ylabel('Total Process Energy (J)');
box on;

plotBarStackGroups(timesStack_norm, groupLabels);
xlabel('Applications');
ylabel('Normalised Process Time');
box on;

yyaxis left;
plotBarStackGroups(energyStack_norm, groupLabels);
ylabel('Normalised Process Energy');
yyaxis right;
plot([1, 2, 3, 4, 5, 6, 7, 8], (timesStack_norm(:, 2) - timesStack_norm(:, 1))*100, '-o')
ylabel('Performance Gain')
xlabel('Applications');
box on;


% totalEnergy_parallel = readEnergy_parallel + writeEnergy_parallel;
% totalEnergy_parallel_base = readEnergy_parallel_base + writeEnergy_parallel_base;
% totalEnergy_serial = readEnergy_serial + writeEnergy_serial;
% totalEnergy_serial_base = readEnergy_serial_base + writeEnergy_serial_base;



x = [1, 2];
%  
% figure;
% hp = bar(x, [mean(times_parallel, 1) * Clk_M_Value, 0]);
% hold on;
% hr = bar(x, [0, mean(times_parallel_base, 1) * Clk_M_Value]);
% errorbar(x, [mean(times_parallel, 1) * Clk_M_Value,  mean(times_parallel_base, 1) * Clk_M_Value], [(mean(times_parallel, 1) - min(times_parallel, [], 1)) * Clk_M_Value, (mean(times_parallel_base, 1) - min(times_parallel_base, [], 1)) * Clk_M_Value], [(max(times_parallel, [], 1) - mean(times_parallel, 1)) * Clk_M_Value, (max(times_parallel_base, [], 1) - mean(times_parallel_base, 1)) * Clk_M_Value], 'Color', 'k', 'LineStyle', 'none');
% hold off;
% xticklabels({'Proposed', 'Reference'});
% ylabel('Write time (seconds)');
% 
% figure;
% hp = bar(x, [mean(readEnergy_parallel, 1), mean(writeEnergy_parallel, 1); 0, 0], 'stacked');
% hold on;
% hr = bar(x, [0, 0; mean(readEnergy_parallel_base, 1), mean(writeEnergy_parallel_base, 1)], 'stacked');
% errorbar(x, [mean(totalEnergy_parallel, 1),  mean(totalEnergy_parallel_base, 1)], [(mean(totalEnergy_parallel, 1) - min(totalEnergy_parallel, [], 1)), (mean(totalEnergy_parallel_base, 1) - min(totalEnergy_parallel_base, [], 1))], [(max(totalEnergy_parallel, [], 1) - mean(totalEnergy_parallel, 1)), (max(totalEnergy_parallel_base, [], 1) - mean(totalEnergy_parallel_base, 1))], 'Color', 'k', 'LineStyle', 'none');
% hold off;
% xticklabels({'Proposed', 'Reference'});
% ylabel('Total Energy (J)');
% 
% figure;
% hp = bar(x, [mean(readStallTimes_parallel, 1), 0]);
% hold on;
% hr = bar(x, [0, mean(readStallTimes_parallel_base, 1)]);
% errorbar(x, [mean(readStallTimes_parallel, 1),  mean(readStallTimes_parallel_base, 1)], [(mean(readStallTimes_parallel, 1) - min(readStallTimes_parallel, [], 1)), (mean(readStallTimes_parallel_base, 1) - min(readStallTimes_parallel_base, [], 1))], [(max(readStallTimes_parallel, [], 1) - mean(readStallTimes_parallel, 1)), (max(readStallTimes_parallel_base, [], 1) - mean(readStallTimes_parallel_base, 1))], 'Color', 'k', 'LineStyle', 'none');
% hold off;
% xticklabels({'Proposed', 'Reference'});
% ylabel('Read Stall Cycles');
% 
% figure;
% hp = bar(x, [mean(writeStallTimes_parallel, 1), 0]);
% hold on;
% hr = bar(x, [0, mean(writeStallTimes_parallel_base, 1)]);
% errorbar(x, [mean(writeStallTimes_parallel, 1),  mean(writeStallTimes_parallel_base, 1)], [(mean(writeStallTimes_parallel, 1) - min(writeStallTimes_parallel, [], 1)), (mean(writeStallTimes_parallel_base, 1) - min(writeStallTimes_parallel_base, [], 1))], [(max(writeStallTimes_parallel, [], 1) - mean(writeStallTimes_parallel, 1)), (max(writeStallTimes_parallel_base, [], 1) - mean(writeStallTimes_parallel_base, 1))], 'Color', 'k', 'LineStyle', 'none');
% hold off;
% xticklabels({'Proposed', 'Reference'});
% ylabel('Write Stall Cycles');


% hold on;
% errorbar(writeBufferSize, mean(times_parallel_base, 1) * Clk_M_Value, (mean(times_parallel_base, 1) - min(times_parallel_base, [], 1)) * Clk_M_Value, (max(times_parallel_base, [], 1) - mean(times_parallel_base, 1)) * Clk_M_Value, '-o');
% hold off;
% legend({sprintf('Parallel Write Proposed (%d bits)', numBits), sprintf('Parallel Write Reference (%d bits)', numBits)})
% savefig('writeTime_vs_batchSize_Parallel.fig');

% figure;
% errorbar(writeBufferSize, mean(readEnergy_parallel, 1), (mean(readEnergy_parallel, 1) - min(readEnergy_parallel, [], 1)), (max(readEnergy_parallel, [], 1) - mean(readEnergy_parallel, 1)), '-o');
% xlabel('Batch Size');
% ylabel('Read Energy (J)');
% 
% hold on;
% errorbar(writeBufferSize, mean(readEnergy_parallel_base, 1), (mean(readEnergy_parallel_base, 1) - min(readEnergy_parallel_base, [], 1)), (max(readEnergy_parallel_base, [], 1) - mean(readEnergy_parallel_base, 1)), '-o');
% hold off;
% legend({sprintf('Parallel Write Proposed (%d bits)', numBits), sprintf('Parallel Write Reference (%d bits)', numBits)})
% % savefig('readEnergy_vs_batchSize_Parallel.fig');
% 
% figure;
% errorbar(writeBufferSize, mean(writeEnergy_parallel, 1), (mean(writeEnergy_parallel, 1) - min(writeEnergy_parallel, [], 1)), (max(writeEnergy_parallel, [], 1) - mean(writeEnergy_parallel, 1)), '-o');
% xlabel('Batch Size');
% ylabel('Write Energy (J)');
% 
% hold on;
% errorbar(writeBufferSize, mean(writeEnergy_parallel_base, 1), (mean(writeEnergy_parallel_base, 1) - min(writeEnergy_parallel_base, [], 1)), (max(writeEnergy_parallel_base, [], 1) - mean(writeEnergy_parallel_base, 1)), '-o');
% hold off;
% legend({sprintf('Parallel Write Proposed (%d bits)', numBits), sprintf('Parallel Write Reference (%d bits)', numBits)})
% % savefig('writeEnergy_vs_batchSize_Parallel.fig');
% 
% figure;
% errorbar(writeBufferSize, mean(totalEnergy_parallel, 1), (mean(totalEnergy_parallel, 1) - min(totalEnergy_parallel, [], 1)), (max(totalEnergy_parallel, [], 1) - mean(totalEnergy_parallel, 1)), '-o');
% xlabel('Batch Size');
% ylabel('Total Energy (J)');
% 
% hold on;
% errorbar(writeBufferSize, mean(totalEnergy_parallel_base, 1), (mean(totalEnergy_parallel_base, 1) - min(totalEnergy_parallel_base, [], 1)), (max(totalEnergy_parallel_base, [], 1) - mean(totalEnergy_parallel_base, 1)), '-o');
% hold off;
% legend({sprintf('Parallel Write Proposed (%d bits)', numBits), sprintf('Parallel Write Reference (%d bits)', numBits)})
% % savefig('totalEnergy_vs_batchSize_Parallel.fig');