% This test gets the stall time contour plots for either a particular 
% application or random inputs.

% You can set different distribution parameters (for e.g. you can change 
% the sigma) or you can change the processor frequency and put multiple
% plots on the same figure. Then, for a prticular set of parameters, you
% can determine the optimal set of waitBuffer and batch size and run the
% performance tests on those size values

clc;
clear all;
close all;

warning('off', 'all');

Clk_M = 1;
Clk_M_Value = 5e-9;
mu = 5; % [5, 5, 5.05];
sigma = 1; % [1, 2, 4]; 
leftDiscard = 0;
rightDiscard = 4; % [4, 2, 1];
numIter = 20;% 200;
numAddresses = 1000;
numBits = 8;
bitReadTime = 1* Clk_M;
switchDetectTime = 0.2 * Clk_M;
global rwArray;
global rwArray_local;

% Offset Multiplier Tests

writeBufferSize = [2, 5:5:100]; % [1:1:4, 5:5:100]; % [2, 5:5:100]; % 1:10;
waitBufferSize = [5:5:120]; % [5:5:100];
[writeBufferSize_fine, waitBufferSize_fine] = meshgrid([2:2:100], [5:2:120]);
offsetMult = [1, 1.05, 1.1];
offsetMultBase = 2; % (mu + rightDiscard * sigma) / mu;
% times_parallel = zeros(size(offsetMult, 2), numIter);
% times_parallel_local = zeros(size(offsetMult, 2), numIter);
% times_parallel_bounded = zeros(size(offsetMult, 2), numIter);
% times_parallel_base = zeros(size(offsetMult, 2), numIter);
% times_serial = zeros(size(offsetMult, 2), numIter);
% times_serial_base = zeros(size(offsetMult, 2), numIter);
% maxBufferSize = zeros(size(offsetMult, 2), numIter);
% maxBufferSize_local = zeros(size(offsetMult, 2), numIter);
stall_times = NaN(length(waitBufferSize), length(writeBufferSize), numIter);
total_times = NaN(length(waitBufferSize), length(writeBufferSize), numIter);

legendInfo = cell(length(waitBufferSize), 1);

tic;
for o = 1:1:length(offsetMult)
    for i = 1:length(waitBufferSize)
        for j = 1:length(writeBufferSize)
            parfor n = 1:numIter
                rwArray = cell(numAddresses, 2);
        %         memTrace = fopen(sprintf('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 6 MMU CS\\ECG_3s_each\\report_global_%d.csv', randi([0, 19])));
        %         % memTrace = fopen(sprintf('/scrap/users/tuli/ECG_3s_each/report_global_%d.csv', randi([0, 19])));
        %         rwArray = {};
        %         % Forming rwArray from memTrace
        %         cellArray = textscan(memTrace, '%c %d %s');
        %         rwArray = cell(length(cellArray{1, 1}), 2);
        %         for k = 1:length(cellArray{1, 1})
        %             rwArray{k, 1} = cellArray{1, 1}(k);
        %             if cellArray{1, 1}(k) == 'w'
        %                 rwArray{k, 2} = [cellArray{1, 2}(k), hex2dec(cellArray{1, 3}(k))];
        %             else
        %                 rwArray{k, 2} = cellArray{1, 2}(k);
        %             end
        %         end
                rwArray_local = {};
    %             for k = 1:1:numAddresses
    %                 rwArray{k, 1} = 'w';
    %                 rwArray{k, 2} = [k, randi([0, 255])];
    %             end
                for k = randi([1, 1 * numAddresses], 1, numAddresses)
                    if size(rwArray_local, 1) > 0
                        rwArray_local{size(rwArray_local, 1) + 1, 1} = 'w';
                        rwArray_local{size(rwArray_local, 1), 2} = [k, randi([0, 255])];
                    else
                        rwArray_local = {'w', [k, randi([0, 255])]};
                    end
                end
                if waitBufferSize(i) >= writeBufferSize(j)
                    output = MMU_rw_parallel_new(rwArray_local, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult(o), writeBufferSize(j), bitReadTime, switchDetectTime, waitBufferSize(i), 0);
                    % output = MMU_rw_parallel(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMultBase, 1, bitReadTime, switchDetectTime, 1, 1);
                    stall_times(i, j, n) = output(1);
                    total_times(i, j, n) = output(2);
                end
        %         times_parallel_bounded(i, j) = MMU_rw_parallel(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult(i), writeBufferSize, bitReadTime, switchDetectTime, waitBufferSize(m), 0);
        %         times_parallel_local(i, j) = MMU_rw_parallel(rwArray_local, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult(i), writeBufferSize, bitReadTime, switchDetectTime, 0);
        %         times_parallel(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult(i), writeBufferSize, 0);
        %         times_serial(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult(i)*8, writeBufferSize, 0);
        %         times_parallel_base(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase, writeBufferSize, 1);
        %         times_serial_base(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase*8, writeBufferSize, 1);
            end
        end
    %     errorbar(2./offsetMult, mean(times_parallel, 2), mean(times_parallel, 2) - min(times_parallel, [], 2), max(times_parallel, [], 2) - mean(times_parallel, 2), '-o');
    %     hold on;
    %     if m == 1
    %         legendInfo{m} = 'Unbounded';
    %     else
    %         legendInfo{m} = sprintf('Maximum waitBuffer size = %d', waitBufferSize(m));
    %     end
    %     legendInfo{m} = sprintf('Maximum waitBuffer size = %d', waitBufferSize(m));
        fprintf('%0.1f%% Completed\n', i/length(waitBufferSize)*100);
    end
    
    stall_times_fine = interp2(writeBufferSize, waitBufferSize, mean(stall_times, 3), writeBufferSize_fine, waitBufferSize_fine, 'linear');
    % contour(writeBufferSize, waitBufferSize, mean(stall_times, 3), [0, 1], 'ShowText', 'off', 'LineColor', 'k');
    contour(writeBufferSize_fine, waitBufferSize_fine, stall_times_fine, [0, 1], 'ShowText', 'off', 'LineColor', 'k');
    xlabel('writeBuffer size');
    ylabel('waitBuffer size');
    title('# Stall cycles');
    hold on;
end
hold off;
toc;

% figure;
% contour(writeBufferSize, waitBufferSize, mean(stall_times, 3), [0, 1], 'ShowText', 'off');
% xlabel('writeBuffer size');
% ylabel('waitBuffer size');
% title('# Stall cycles');
% hold on;

% figure;
% heatmap(writeBufferSize, flip(waitBufferSize, 2), flip(mean(total_times, 3), 1));
% xlabel('writeBuffer size');
% ylabel('waitBuffer size');
% title('# Total time (Clk_Ms)');

% figure;
% errorbar(2./offsetMult, mean(times_parallel, 2), mean(times_parallel, 2) - min(times_parallel, [], 2), max(times_parallel, [], 2) - mean(times_parallel, 2), '-o');
% % xlabel('Offset Multiplier');
% xlabel('f_P / f_{M(WC)}');
% ylabel('Write time (Clk_Ms)');
% hold on;
% errorbar(2./offsetMult, mean(times_parallel_bounded, 2), mean(times_parallel_bounded, 2) - min(times_parallel_bounded, [], 2), max(times_parallel_bounded, [], 2) - mean(times_parallel_bounded, 2), '-o');
% legend({'Unbounded waitBuffer', 'Bounded waitBuffer'});

