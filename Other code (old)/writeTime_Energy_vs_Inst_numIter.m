clc;
close all;
clear all;


mu = 5;
sigma = 2; 
leftDiscard = 0;
rightDiscard = 2;
numBits = 32;
Clk_M = 1;
Clk_M_Value = 5e-9; % in seconds
bitReadTime = 1* Clk_M;
switchDetectTime = 0.2 * Clk_M; % 1ns
writeVolt = 1; % in Volts
I_LRS = 100e-6; % in Amps
I_HRS = 15e-6; % in Amps
readEnergyPerBit = 1e-12; % in J
debug = 0;
rwArray = {};

writeBufferSize = 10;
numIter = 5;
offsetMult = 1.1; % 1:0.1:2;
offsetMultBase = (mu + rightDiscard * sigma) / mu;
times_parallel = zeros(size(writeBufferSize, 2), numIter);
times_parallel_base = zeros(size(writeBufferSize, 2), numIter);
times_serial = zeros(size(writeBufferSize, 2), numIter);
times_serial_base = zeros(size(writeBufferSize, 2), numIter);
readEnergy_parallel = zeros(size(writeBufferSize, 2), numIter);
readEnergy_parallel_base = zeros(size(writeBufferSize, 2), numIter);
readEnergy_serial = zeros(size(writeBufferSize, 2), numIter);
readEnergy_serial_base = zeros(size(writeBufferSize, 2), numIter); 
writeEnergy_parallel = zeros(size(writeBufferSize, 2), numIter);
writeEnergy_parallel_base = zeros(size(writeBufferSize, 2), numIter);
writeEnergy_serial = zeros(size(writeBufferSize, 2), numIter);
writeEnergy_serial_base = zeros(size(writeBufferSize, 2), numIter);

tic;
for i = 1:length(writeBufferSize)
    for j = 1:numIter
        if i == 1 && j == 1
            memTrace = fopen(sprintf('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 6 MMU CS\\ECG_3s_each\\report_global_%d.csv', randi([0, 19])));
            % memTrace = fopen(sprintf('/scrap/users/tuli/ECG_3s_each/report_global_%d.csv', randi([0, 19])));
            rwArray = {};
            % Forming rwArray from memTrace
            cellArray = textscan(memTrace, '%c %d %s');
            rwArray = cell(5000, 2); % cell(length(cellArray{1, 1}), 2);
            for k = 1:5000 % length(cellArray{1, 1})
                rwArray{k, 1} = cellArray{1, 1}(k);
                if cellArray{1, 1}(k) == 'w'
                    rwArray{k, 2} = [cellArray{1, 2}(k), hex2dec(cellArray{1, 3}(k))];
                else
                    rwArray{k, 2} = cellArray{1, 2}(k);
                end
            end
        end
        writeTimes_inst = zeros(size(rwArray, 1), numIter);
        writeEnergy_inst = zeros(size(rwArray, 1), numIter);
        writeTimes_inst_base = zeros(size(rwArray, 1), numIter);
        writeEnergy_inst_base = zeros(size(rwArray, 1), numIter);
        [times_parallel(i, j), readEnergy_parallel(i, j), writeEnergy_parallel(i, j), writeTimes_inst(:, j), writeEnergy_inst(:, j)] = MMU_rw_parallel_wEnergy_wInst(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult, writeBufferSize(i), bitReadTime, switchDetectTime, writeVolt, I_LRS, I_HRS, readEnergyPerBit, Clk_M_Value, 0, debug);
        % [times_serial(i, j), readEnergy_serial(i, j), writeEnergy_serial(i, j)] = MMU_rw_serial_wEnergy(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult*8, writeBufferSize(i), bitReadTime, switchDetectTime, writeVolt, I_LRS, I_HRS, readEnergyPerBit, Clk_M_Value, 0, debug);
        [times_parallel_base(i, j), readEnergy_parallel_base(i, j), writeEnergy_parallel_base(i, j), writeTimes_inst_base(:, j), writeEnergy_inst_base(:, j)] = MMU_rw_parallel_wEnergy_wInst(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMultBase, writeBufferSize(i), bitReadTime, switchDetectTime, writeVolt, I_LRS, I_HRS, readEnergyPerBit, Clk_M_Value, 1, debug);
        % [times_serial_base(i, j), readEnergy_serial_base(i, j), writeEnergy_serial_base(i, j)] = MMU_rw_serial_wEnergy(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMultBase*8, writeBufferSize(i), bitReadTime, switchDetectTime, writeVolt, I_LRS, I_HRS, readEnergyPerBit, Clk_M_Value, 1, debug);
    end
    fprintf('%0.1f%% Completed\n', i/length(writeBufferSize)*100);
    toc;
end

% writeTimes_inst_new = [];
% inst1 = [];
% inst1_idx = 0;
% writeTimes_inst_base_new = [];
% inst2 = [];
% inst2_isd = 0;
% writeEnergy_inst_new = [];
% inst3 = [];
% inst3_idx = 0;
% writeEnergy_inst_base_new = [];
% inst4 = [];
% inst4_idx = 0;
% 
% for i = 1:1:size(rwArray, 1)
%     if writeTimes_inst(i) ~= 0
%         writeTimes_inst_new = [writeTimes_inst_new; writeTimes_inst(i)];
%         inst1 = [inst1; i];
%     end
% end
% 
% for i = 1:1:size(rwArray, 1)
%     if writeTimes_inst_base(i) ~= 0
%         writeTimes_inst_base_new = [writeTimes_inst_base_new; writeTimes_inst_base(i)];
%         inst2 = [inst2; i];
%     end
% end
% 
% for i = 1:1:size(rwArray, 1)
%     if writeEnergy_inst(i) ~= 0
%         writeEnergy_inst_new = [writeEnergy_inst_new; writeEnergy_inst(i)];
%         inst3 = [inst3; i];
%     end
% end
% 
% for i = 1:1:size(rwArray, 1)
%     if writeEnergy_inst_base(i) ~= 0
%         writeEnergy_inst_base_new = [writeEnergy_inst_base_new; writeEnergy_inst_base(i)];
%         inst4 = [inst4; i];
%     end
% end


figure;
errorbar(1:1:size(rwArray, 1), mean(writeTimes_inst(1:1:size(rwArray, 1), :), 2) * Clk_M_Value, (mean(writeTimes_inst(1:1:size(rwArray, 1), :), 2) - min(writeTimes_inst(1:1:size(rwArray, 1), :), [], 2)) * Clk_M_Value, (max(writeTimes_inst(1:1:size(rwArray, 1), :), [], 2) - mean(writeTimes_inst(1:1:size(rwArray, 1), :), 2)) * Clk_M_Value, 'o')
% scatter(inst1(51:end), writeTimes_inst_new(51:end) * Clk_M_Value, 5, 'filled')
hold on
% plot(inst1(51:end), movmean(writeTimes_inst_new(51:end) * Clk_M_Value, [5, 0]), '--k')
errorbar(1:1:size(rwArray, 1), mean(writeTimes_inst_base(1:1:size(rwArray, 1), :), 2) * Clk_M_Value, (mean(writeTimes_inst_base(1:1:size(rwArray, 1), :), 2) - min(writeTimes_inst_base(1:1:size(rwArray, 1), :), [], 2)) * Clk_M_Value, (max(writeTimes_inst_base(1:1:size(rwArray, 1), :), [], 2) - mean(writeTimes_inst_base(1:1:size(rwArray, 1), :), 2)) * Clk_M_Value, 'o')
% scatter(inst2(51:end), writeTimes_inst_base_new(51:end) * Clk_M_Value, 5, 'filled')
% plot(inst2(51:end), movmean(writeTimes_inst_base_new(51:end) * Clk_M_Value, 5), '--r')
hold off
xlabel('Instructions')
ylabel('Write Times (seconds)')
legend({sprintf('Parallel Write Proposed (%d bits)', numBits), sprintf('Parallel Write Reference (%d bits)', numBits)});

% figure;
% scatter(1:1:size(rwArray, 1), writeEnergy_inst, '.')
% % scatter(inst3(51:end), writeEnergy_inst_new(51:end), 5, 'filled')
% hold on
% % plot(inst3(51:end), movmean(writeEnergy_inst_new(51:end), [5, 0]), '--k')
% scatter(1:1:size(rwArray, 1), writeEnergy_inst_base, '.')
% % scatter(inst4(51:end), writeEnergy_inst_base_new(51:end), 5, 'filled')
% % plot(inst4(51:end), movmean(writeEnergy_inst_base_new(51:end), [5, 0]), '--k')
% hold off
% xlabel('Instructions')
% ylabel('Write Energy (J)')
% legend({sprintf('Parallel Write Proposed (%d bits)', numBits), 'Moving Average', sprintf('Parallel Write Reference (%d bits)', numBits), 'Moving Average'});

% totalEnergy_parallel = readEnergy_parallel + writeEnergy_parallel;
% totalEnergy_parallel_base = readEnergy_parallel_base + writeEnergy_parallel_base;
% totalEnergy_serial = readEnergy_serial + writeEnergy_serial;
% totalEnergy_serial_base = readEnergy_serial_base + writeEnergy_serial_base;
% 
% figure;
% errorbar(writeBufferSize, mean(times_parallel, 2) * Clk_M_Value, (mean(times_parallel, 2) - min(times_parallel, [], 2)) * Clk_M_Value, (max(times_parallel, [], 2) - mean(times_parallel, 2)) * Clk_M_Value, '-o');
% xlabel('Batch Size');
% ylabel('Write time (seconds)');
% 
% hold on;
% errorbar(writeBufferSize, mean(times_parallel_base, 2) * Clk_M_Value, (mean(times_parallel_base, 2) - min(times_parallel_base, [], 2)) * Clk_M_Value, (max(times_parallel_base, [], 2) - mean(times_parallel_base, 2)) * Clk_M_Value, '-o');
% hold off;
% legend({sprintf('Parallel Write Proposed (%d bits)', numBits), sprintf('Parallel Write Reference (%d bits)', numBits)})
% savefig('writeTime_vs_batchSize_Parallel.fig');
% 
% figure;
% errorbar(writeBufferSize, mean(times_serial, 2) * Clk_M_Value, (mean(times_serial, 2) - min(times_serial, [], 2)) * Clk_M_Value, (max(times_serial, [], 2) - mean(times_serial, 2)) * Clk_M_Value, '-o');
% xlabel('Batch Size');
% ylabel('Write time (seconds)');
% 
% hold on;
% errorbar(writeBufferSize, mean(times_serial_base, 2) * Clk_M_Value, (mean(times_serial_base, 2) - min(times_serial_base, [], 2)) * Clk_M_Value, (max(times_serial_base, [], 2) - mean(times_serial_base, 2)) * Clk_M_Value, '-o');
% hold off;
% legend({'Serial Write Proposed', 'Serial Write Reference'})
% savefig('writeTime_vs_batchSize_Serial.fig');
% 
% figure;
% errorbar(writeBufferSize, mean(readEnergy_parallel, 2), (mean(readEnergy_parallel, 2) - min(readEnergy_parallel, [], 2)), (max(readEnergy_parallel, [], 2) - mean(readEnergy_parallel, 2)), '-o');
% xlabel('Batch Size');
% ylabel('Read Energy (J)');
% 
% hold on;
% errorbar(writeBufferSize, mean(readEnergy_parallel_base, 2), (mean(readEnergy_parallel_base, 2) - min(readEnergy_parallel_base, [], 2)), (max(readEnergy_parallel_base, [], 2) - mean(readEnergy_parallel_base, 2)), '-o');
% hold off;
% legend({sprintf('Parallel Write Proposed (%d bits)', numBits), sprintf('Parallel Write Reference (%d bits)', numBits)})
% savefig('readEnergy_vs_batchSize_Parallel.fig');
% 
% figure;
% errorbar(writeBufferSize, mean(readEnergy_serial, 2), (mean(readEnergy_serial, 2) - min(readEnergy_serial, [], 2)), (max(readEnergy_serial, [], 2) - mean(readEnergy_serial, 2)), '-o');
% xlabel('Batch Size');
% ylabel('Read Energy (J)');
% 
% hold on;
% errorbar(writeBufferSize, mean(readEnergy_serial_base, 2), (mean(readEnergy_serial_base, 2) - min(readEnergy_serial_base, [], 2)), (max(readEnergy_serial_base, [], 2) - mean(readEnergy_serial_base, 2)), '-o');
% hold off;
% legend({'Serial Write Proposed', 'Serial Write Reference'})
% savefig('readEnergy_vs_batchSize_Serial.fig');
% 
% figure;
% errorbar(writeBufferSize, mean(writeEnergy_parallel, 2), (mean(writeEnergy_parallel, 2) - min(writeEnergy_parallel, [], 2)), (max(writeEnergy_parallel, [], 2) - mean(writeEnergy_parallel, 2)), '-o');
% xlabel('Batch Size');
% ylabel('Write Energy (J)');
% 
% hold on;
% errorbar(writeBufferSize, mean(writeEnergy_parallel_base, 2), (mean(writeEnergy_parallel_base, 2) - min(writeEnergy_parallel_base, [], 2)), (max(writeEnergy_parallel_base, [], 2) - mean(writeEnergy_parallel_base, 2)), '-o');
% hold off;
% legend({sprintf('Parallel Write Proposed (%d bits)', numBits), sprintf('Parallel Write Reference (%d bits)', numBits)})
% savefig('writeEnergy_vs_batchSize_Parallel.fig');
% 
% figure;
% errorbar(writeBufferSize, mean(writeEnergy_serial, 2), (mean(writeEnergy_serial, 2) - min(writeEnergy_serial, [], 2)), (max(writeEnergy_serial, [], 2) - mean(writeEnergy_serial, 2)), '-o');
% xlabel('Batch Size');
% ylabel('Write Energy (J)');
% 
% hold on;
% errorbar(writeBufferSize, mean(writeEnergy_serial_base, 2), (mean(writeEnergy_serial_base, 2) - min(writeEnergy_serial_base, [], 2)), (max(writeEnergy_serial_base, [], 2) - mean(writeEnergy_serial_base, 2)), '-o');
% hold off;
% legend({'Serial Write Proposed', 'Serial Write Reference'})
% savefig('writeEnergy_vs_batchSize_Serial.fig');
% 
% figure;
% errorbar(writeBufferSize, mean(totalEnergy_parallel, 2), (mean(totalEnergy_parallel, 2) - min(totalEnergy_parallel, [], 2)), (max(totalEnergy_parallel, [], 2) - mean(totalEnergy_parallel, 2)), '-o');
% xlabel('Batch Size');
% ylabel('Total Energy (J)');
% 
% hold on;
% errorbar(writeBufferSize, mean(totalEnergy_parallel_base, 2), (mean(totalEnergy_parallel_base, 2) - min(totalEnergy_parallel_base, [], 2)), (max(totalEnergy_parallel_base, [], 2) - mean(totalEnergy_parallel_base, 2)), '-o');
% hold off;
% legend({sprintf('Parallel Write Proposed (%d bits)', numBits), sprintf('Parallel Write Reference (%d bits)', numBits)})
% savefig('totalEnergy_vs_batchSize_Parallel.fig');
% 
% figure;
% errorbar(writeBufferSize, mean(totalEnergy_serial, 2), (mean(totalEnergy_serial, 2) - min(totalEnergy_serial, [], 2)), (max(totalEnergy_serial, [], 2) - mean(totalEnergy_serial, 2)), '-o');
% xlabel('Batch Size');
% ylabel('Total Energy (J)');
% 
% hold on;
% errorbar(writeBufferSize, mean(totalEnergy_serial_base, 2), (mean(totalEnergy_serial_base, 2) - min(totalEnergy_serial_base, [], 2)), (max(totalEnergy_serial_base, [], 2) - mean(totalEnergy_serial_base, 2)), '-o');
% hold off;
% legend({'Serial Write Proposed', 'Serial Write Reference'})
% savefig('totalEnergy_vs_batchSize_Serial.fig');