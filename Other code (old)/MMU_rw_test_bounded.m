clc;
clear all;
% close all;

warning('off', 'all');

Clk_M = 1;
Clk_M_Value = 5e-9;
mu = 5;
sigma = 2; 
leftDiscard = 0;
rightDiscard = 2;
numIter = 100;
numAddresses = 100;
numBits = 32;
bitReadTime = 1* Clk_M;
switchDetectTime = 0.2 * Clk_M;
global rwArray;

% Offset Multiplier Tests

writeBufferSize = 5; % 1:10;
maxWaitBufferSize = [10:10:50];
offsetMult = 1:0.1:2; % 1.1;
offsetMultBase = 2; % (mu + rightDiscard * sigma) / mu;
times_parallel = zeros(size(offsetMult, 2), numIter);
% times_parallel_local = zeros(size(offsetMult, 2), numIter);
times_parallel_bounded = zeros(size(offsetMult, 2), numIter);
times_parallel_base = zeros(size(offsetMult, 2), numIter);
times_serial = zeros(size(offsetMult, 2), numIter);
times_serial_base = zeros(size(offsetMult, 2), numIter);
% maxBufferSize = zeros(size(offsetMult, 2), numIter);
% maxBufferSize_local = zeros(size(offsetMult, 2), numIter);

legendInfo = cell(length(maxWaitBufferSize), 1);

figure;

for m = 1:length(maxWaitBufferSize)
    for i = 1:length(offsetMult)
        parfor j = 1:numIter
            rwArray = {};
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
    %         rwArray_local = {};
            for k = randperm(numAddresses)
                if size(rwArray, 1) > 0
                    rwArray{size(rwArray, 1) + 1, 1} = 'w';
                    rwArray{size(rwArray, 1), 2} = [k, randi([0, 255])];
                else
                    rwArray = {'w', [k, randi([0, 255])]};
                end
            end
    %         for k = randi([1, 0.75 * numAddresses], 1, numAddresses)
    %             if size(rwArray_local, 1) > 0
    %                 rwArray_local{size(rwArray_local, 1) + 1, 1} = 'w';
    %                 rwArray_local{size(rwArray_local, 1), 2} = [k, randi([0, 255])];
    %             else
    %                 rwArray_local = {'w', [k, randi([0, 255])]};
    %             end
    %         end
            times_parallel(i, j) = MMU_rw_parallel(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult(i), writeBufferSize, bitReadTime, switchDetectTime, maxWaitBufferSize(m), 0);
    %         times_parallel_bounded(i, j) = MMU_rw_parallel(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult(i), writeBufferSize, bitReadTime, switchDetectTime, maxWaitBufferSize(m), 0);
    %         times_parallel_local(i, j) = MMU_rw_parallel(rwArray_local, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult(i), writeBufferSize, bitReadTime, switchDetectTime, 0);
    %         times_parallel(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult(i), writeBufferSize, 0);
    %         times_serial(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult(i)*8, writeBufferSize, 0);
    %         times_parallel_base(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase, writeBufferSize, 1);
    %         times_serial_base(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase*8, writeBufferSize, 1);
        end
    end
    errorbar(2./offsetMult, mean(times_parallel, 2), mean(times_parallel, 2) - min(times_parallel, [], 2), max(times_parallel, [], 2) - mean(times_parallel, 2), '-o');
    hold on;
%     if m == 1
%         legendInfo{m} = 'Unbounded';
%     else
%         legendInfo{m} = sprintf('Maximum waitBuffer size = %d', maxWaitBufferSize(m));
%     end
    legendInfo{m} = sprintf('Maximum waitBuffer size = %d', maxWaitBufferSize(m));
    fprintf('%0.1f%% Completed\n', m/length(maxWaitBufferSize)*100);
end

hold off;
xlabel('f_P / f_{M(WC)}');
ylabel('Write time (Clk_Ms)');
legend(legendInfo);

% figure;
% errorbar(2./offsetMult, mean(times_parallel, 2), mean(times_parallel, 2) - min(times_parallel, [], 2), max(times_parallel, [], 2) - mean(times_parallel, 2), '-o');
% % xlabel('Offset Multiplier');
% xlabel('f_P / f_{M(WC)}');
% ylabel('Write time (Clk_Ms)');
% hold on;
% errorbar(2./offsetMult, mean(times_parallel_bounded, 2), mean(times_parallel_bounded, 2) - min(times_parallel_bounded, [], 2), max(times_parallel_bounded, [], 2) - mean(times_parallel_bounded, 2), '-o');
% legend({'Unbounded waitBuffer', 'Bounded waitBuffer'});

