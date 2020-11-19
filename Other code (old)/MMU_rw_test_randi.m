clc;
clear all;
% close all;

warning('off', 'all');

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

numIter = 30;
numAddresses = 20000;
numBits = 8;
bitReadTime = 1* Clk_M;
switchDetectTime = 0.2 * Clk_M;
global rwArray;

% Offset Multiplier Tests

% writeBufferSize = 5; % 1:10;
% offsetMult = 1:0.1:2; % 1.1;
% offsetMultBase = 2; % (mu + rightDiscard * sigma) / mu;
% times_parallel = zeros(size(offsetMult, 2), numIter);
% times_parallel_local = zeros(size(offsetMult, 2), numIter);
% times_parallel_bounded = zeros(size(offsetMult, 2), numIter);
% times_parallel_base = zeros(size(offsetMult, 2), numIter);
% times_serial = zeros(size(offsetMult, 2), numIter);
% times_serial_base = zeros(size(offsetMult, 2), numIter);
% maxBufferSize = zeros(size(offsetMult, 2), numIter);
% maxBufferSize_local = zeros(size(offsetMult, 2), numIter);

% localAdd = [];
% for i = 1:1:0.5*numAddresses
%     localAdd = [localAdd, i, i];
% end

% for s = 1:length(sigma)
% for i = 1:length(offsetMult)
%     parfor j = 1:numIter
%         rwArray = {};
% %         memTrace = fopen(sprintf('C:\\Users\\tuli\\Desktop\\Shikhar Tuli ESL\\Results 6 MMU CS\\ECG_3s_each\\report_global_%d.csv', randi([0, 19])));
% %         % memTrace = fopen(sprintf('/scrap/users/tuli/ECG_3s_each/report_global_%d.csv', randi([0, 19])));
% %         rwArray = {};
% %         % Forming rwArray from memTrace
% %         cellArray = textscan(memTrace, '%c %d %s');
% %         rwArray = cell(length(cellArray{1, 1}), 2);
% %         for k = 1:length(cellArray{1, 1})
% %             rwArray{k, 1} = cellArray{1, 1}(k);
% %             if cellArray{1, 1}(k) == 'w'
% %                 rwArray{k, 2} = [cellArray{1, 2}(k), hex2dec(cellArray{1, 3}(k))];
% %             else
% %                 rwArray{k, 2} = cellArray{1, 2}(k);
% %             end
% %         end
% %         rwArray_local = {};
%         for k = randperm(numAddresses)
%             if size(rwArray, 1) > 0
%                 rwArray{size(rwArray, 1) + 1, 1} = 'w';
%                 rwArray{size(rwArray, 1), 2} = [k, randi([0, 255])];
%             else
%                 rwArray = {'w', [k, randi([0, 255])]};
%             end
%         end
% %         for k = randi([1, 0.75 * numAddresses], 1, numAddresses)
% %             if size(rwArray_local, 1) > 0
% %                 rwArray_local{size(rwArray_local, 1) + 1, 1} = 'w';
% %                 rwArray_local{size(rwArray_local, 1), 2} = [k, randi([0, 255])];
% %             else
% %                 rwArray_local = {'w', [k, randi([0, 255])]};
% %             end
% %         end
%         output = MMU_rw_parallel(rwArray, mu, sigma(s), numBits, leftDiscard, rightDiscard, offsetMult(i), writeBufferSize, bitReadTime, switchDetectTime, 0, 0);
%         times_parallel(i, j) = output(3);
% %         times_parallel_bounded(i, j) = MMU_rw_parallel(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult(i), writeBufferSize, bitReadTime, switchDetectTime, 1, 0);
% %         times_parallel_local(i, j) = MMU_rw_parallel(rwArray_local, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult(i), writeBufferSize, bitReadTime, switchDetectTime, 0);
% %         times_parallel(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult(i), writeBufferSize, 0);
% %         times_serial(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult(i)*8, writeBufferSize, 0);
% %         times_parallel_base(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase, writeBufferSize, 1);
% %         times_serial_base(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase*8, writeBufferSize, 1);
%     end
%     fprintf('%0.1f%% Completed\n', i/length(offsetMult)*100);
% end
% errorbar(2./offsetMult, mean(times_parallel, 2), mean(times_parallel, 2) - min(times_parallel, [], 2), max(times_parallel, [], 2) - mean(times_parallel, 2), '-o');
% xlabel('f_P / f_{M(WC)}');
% ylabel('Maximum Buffer Size');
% hold on;
% end
% legend({'\sigma = 1', '\sigma = 2', '\sigma = 4'});
% hold off;

% figure;
% errorbar(2./offsetMult, mean(times_parallel, 2), mean(times_parallel, 2) - min(times_parallel, [], 2), max(times_parallel, [], 2) - mean(times_parallel, 2), '-o');
% xlabel('Offset Multiplier');
% xlabel('f_P / f_{M(WC)}');
% ylabel('Maximum Buffer Size');
% hold on;
% errorbar(2./offsetMult, mean(times_parallel_bounded, 2), mean(times_parallel_bounded, 2) - min(times_parallel_bounded, [], 2), max(times_parallel_bounded, [], 2) - mean(times_parallel_bounded, 2), '-o');
% legend({'Unbounded waitBuffer', 'Bounded waitBuffer'});

% figure;
% errorbar(offsetMult, mean(times_parallel, 2), mean(times_parallel, 2) - min(times_parallel, [], 2), max(times_parallel, [], 2) - mean(times_parallel, 2), '-o');
% xlabel('Offset Multiplier');
% ylabel('Write time (Clk_Ms)');
% 
% hold on;
% errorbar(offsetMult, mean(times_parallel_base, 2), mean(times_parallel_base, 2) - min(times_parallel_base, [], 2), max(times_parallel_base, [], 2) - mean(times_parallel_base, 2), '-o');
% hold off;
% legend({'Parallel Write Proposed (8 bits)', 'Parallel Write Reference (8 bits)'})
% 
% figure;
% errorbar(offsetMult, mean(times_serial, 2), mean(times_serial, 2) - min(times_serial, [], 2), max(times_serial, [], 2) - mean(times_serial, 2), '-o');
% xlabel('Offset Multiplier');
% ylabel('Write time (Clk_Ms)');
% 
% hold on;
% errorbar(offsetMult, mean(times_serial_base, 2), mean(times_serial_base, 2) - min(times_serial_base, [], 2), max(times_serial_base, [], 2) - mean(times_serial_base, 2), '-o');
% hold off;
% legend({'Serial Write Proposed', 'Serial Write Reference'})


% Buffer Size Test

% rwSize = numAddresses + 1; % rwSize = numAddresses + 1 means that only write operations performed
% writeBufferSize = [1:1:4, 5:5:50]; % 2;
% offsetMult = 1; % 1:0.1:2;
% offsetMultBase = 2; % (mu + rightDiscard * sigma) / mu;
% times_parallel = zeros(size(writeBufferSize, 2), numIter);
% times_parallel_base = zeros(size(writeBufferSize, 2), numIter);
% times_serial = zeros(size(writeBufferSize, 2), numIter);
% times_serial_base = zeros(size(writeBufferSize, 2), numIter);
% times_parallel = zeros(size(writeBufferSize, 2), numIter);
% times_parallel_local = zeros(size(writeBufferSize, 2), numIter);
% 
% for i = 1:length(writeBufferSize)
%     parfor j = 1:numIter
%         rwArray = {};
%         rwArray_local = {};
%         for k = randperm(numAddresses)
%             if size(rwArray, 1) > 0
%                 rwArray{size(rwArray, 1) + 1, 1} = 'w';
%                 rwArray{size(rwArray, 1), 2} = [k, randi([0, 255])];
%             else
%                 rwArray = {'w', [k, randi([0, 255])]};
%             end
%         end
%         for k = randi([1, 0.75 * numAddresses], 1, numAddresses)
%             if size(rwArray_local, 1) > 0
%                 rwArray_local{size(rwArray_local, 1) + 1, 1} = 'w';
%                 rwArray_local{size(rwArray_local, 1), 2} = [k, randi([0, 255])];
%             else
%                 rwArray_local = {'w', [k, randi([0, 255])]};
%             end
%         end
%         
%         times_parallel(i, j) = MMU_rw_parallel(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult, writeBufferSize(i), bitReadTime, switchDetectTime, 0);
%         times_parallel_local(i, j) = MMU_rw_parallel(rwArray_local, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult, writeBufferSize(i), bitReadTime, switchDetectTime, 0);
%         count = 0;
%         flag = 0;
%         w = 1;
%         for k = 1:numAddresses
%             count = count + 1;
%             if count == rwSize
%                 w = 1 - w;
%                 count = 0;
%             end
%             if w == 1
%                 rwArray{size(rwArray, 1) + 1, 1} = 'w';
%                 rwArray{size(rwArray, 1), 2} = [randi([1, numAddresses]), randi([0, 255])];
% %                 flag = 1 - flag;
% %                 rwArray{size(rwArray, 1), 2} = [1, randi([0, 255])];
% %                 rn = randn;
% %                 while rn  < -1 || rn > 1
% %                     rn = randn;
% %                 end
% %                 rwArray{size(rwArray, 1), 2} = [randi([1, 10]), bin2dec(strcat(dec2bin(round(7.5 + 0.5*rn)), dec2bin(randi([0, 15]))))];
% %                 rwArray{size(rwArray, 1), 2} = [1, bin2dec(strcat(dec2bin(round(7.5 + 0.5*rn)), dec2bin(randi([0, 15]))))];
% %                 if flag == 1
% %                     rwArray{size(rwArray, 1), 2} = [1, 0];
% %                 else
% %                     rwArray{size(rwArray, 1), 2} = [1, 255];
% %                 end
%             else
%                 rwArray{size(rwArray, 1) + 1, 1} = 'r';
%                 rwArray{size(rwArray, 1), 2} = randi([1, numAddresses]);
%             end
%         end
%         times_parallel(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult, writeBufferSize(i), 0);
%         times_serial(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult*8, writeBufferSize(i), 0);
%         times_parallel_base(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase, 1, 1);
%         times_serial_base(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase*8, 1, 1);
%     end
%     fprintf('%0.1f%% Completed\n', i/length(writeBufferSize)*100);
% end
% 
% figure;
% errorbar(writeBufferSize, mean(times_parallel, 2) * Clk_M_Value, (mean(times_parallel, 2) - min(times_parallel, [], 2)) * Clk_M_Value, (max(times_parallel, [], 2) - mean(times_parallel, 2)) * Clk_M_Value, '-o');
% xlabel('Batch Size');
% % xlabel('f_P / f_{M(WC)}');
% ylabel('Write Time (seconds)');
% hold on;
% errorbar(writeBufferSize, mean(times_parallel_local, 2) * Clk_M_Value, (mean(times_parallel_local, 2) - min(times_parallel_local, [], 2)) * Clk_M_Value, (max(times_parallel_local, [], 2) - mean(times_parallel_local, 2)) * Clk_M_Value, '-o');
% legend({'Random addressing', '25% Local addressing'});

% figure;
% errorbar(writeBufferSize, mean(times_parallel, 2), mean(times_parallel, 2) - min(times_parallel, [], 2), max(times_parallel, [], 2) - mean(times_parallel, 2), '-o');
% xlabel('Batch Size');
% ylabel('Write time (Clk_Ms)');
% 
% hold on;
% errorbar(writeBufferSize, mean(times_parallel_base, 2), mean(times_parallel_base, 2) - min(times_parallel_base, [], 2), max(times_parallel_base, [], 2) - mean(times_parallel_base, 2), '-o');
% hold off;
% legend({'Parallel Write Proposed (8 bits)', 'Parallel Write Reference (8 bits)'})
% 
% figure;
% errorbar(writeBufferSize, mean(times_serial, 2), mean(times_serial, 2) - min(times_serial, [], 2), max(times_serial, [], 2) - mean(times_serial, 2), '-o');
% xlabel('Batch Size');
% ylabel('Write time (Clk_Ms)');
% 
% hold on;
% errorbar(writeBufferSize, mean(times_serial_base, 2), mean(times_serial_base, 2) - min(times_serial_base, [], 2), max(times_serial_base, [], 2) - mean(times_serial_base, 2), '-o');
% hold off;
% legend({'Serial Write Proposed', 'Serial Write Reference'})

% Corner cases Test
% 
% % rwSize = numAddresses + 1;
% writeBufferSize = 1; % 2;
% offsetMult = 1.1; % 1:0.1:2;
% offsetMultBase = (mu + rightDiscard * sigma) / mu;
% times_parallel = zeros(2, numIter);
% times_parallel_base = zeros(2, numIter);
% times_serial = zeros(2, numIter);
% times_serial_base = zeros(2, numIter);
% 
% for i = 1:2
%     parfor j = 1:numIter
%         rwArray = {};
%         count = 0;
%         flag = 0;
%         w = 1;
%         for k = 1:numAddresses
%             rwArray{size(rwArray, 1) + 1, 1} = 'w';
% %             rwArray{size(rwArray, 1), 2} = [randi([1, numAddresses]), randi([0, 255])];
%             flag = 1 - flag;
%             if flag == 1
%                 rwArray{size(rwArray, 1), 2} = [1, 0];
%             else
%                 if i == 1
%                     rwArray{size(rwArray, 1), 2} = [1, 0];
%                 else
%                     rwArray{size(rwArray, 1), 2} = [1, 255];
%                 end
%             end
%         end
%         times_parallel(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult, i, 0);
%         times_serial(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult*8, i, 0);
%         times_parallel_base(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase, 1, 1);
%         times_serial_base(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase*8, 1, 1);
%     end
%     fprintf('%d%% Completed\n', i/2*100);
% end
% 
% data = [mean(times_parallel, 2), mean(times_parallel_base, 2)];
% 
% figure;
% bar(data);
% hold on;
% 
% errNeg = [mean(times_parallel, 2) - min(times_parallel, [], 2), mean(times_parallel_base, 2) - min(times_parallel_base, [], 2)];
% errPos = [max(times_parallel, [], 2) - mean(times_parallel, 2),  max(times_parallel_base, [], 2) - mean(times_parallel_base)];
% 
% ngroups = 2;
% nbars = 2;
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% 
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, data(:,i), errNeg(:,i), errPos(:, i), 'k', 'linestyle', 'none');
% end
% hold off;
% set(gca, 'xticklabel', categorical({'0-0 Write op', '0-255-0 Write op'}));
% ylabel('writeBufferSize (Clk_Ms)');
% 
% data = [mean(times_serial, 2), mean(times_serial_base, 2)];
% 
% figure;
% bar(data);
% hold on;
% 
% errNeg = [mean(times_serial, 2) - min(times_serial, [], 2), mean(times_serial_base, 2) - min(times_serial_base, [], 2)];
% errPos = [max(times_serial, [], 2) - mean(times_serial, 2),  max(times_serial_base, [], 2) - mean(times_serial_base)];
% 
% ngroups = 2;
% nbars = 2;
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% 
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, data(:,i), errNeg(:,i), errPos(:, i), 'k', 'linestyle', 'none');
% end
% hold off;
% set(gca, 'xticklabel', categorical({'0-0 Write op', '0-255-0 Write op'}));
% ylabel('writeTime (Clk_Ms)');

% rwSize Test

rwSize = [1, 2, 8, 16, 32];
writeBufferSize = 10; % 1:10;
maxWaitBufferSize = 80;
offsetMult = 1; % 1:0.1:2;
offsetMultBase = 2; % (mu + rightDiscard * sigma) / mu;
times_parallel = zeros(size(rwSize, 2), numIter);
times_parallel_base = zeros(size(rwSize, 2), numIter);
numBits = 8;

% times_serial = zeros(size(rwSize, 2), numIter);
% times_serial_base = zeros(size(rwSize, 2), numIter);

for i = 1:length(rwSize)
    for j = 1:numIter
        rwArray = {};
        count = 0;
        w = 1;
        for k = 1:numAddresses
            count = count + 1;
            if count == rwSize(i)
                w = 1 - w;
                count = 0;
            end
            if w == 1
                rwArray{size(rwArray, 1) + 1, 1} = 'w';
                rwArray{size(rwArray, 1), 2} = [randi([1, numAddresses]), randi([0, 255])];
%                 rwArray{size(rwArray, 1), 2} = [1, randi([0, 255])];
%                 rn = randn;
%                 while rn  < -1 || rn > 1
%                     rn = randn;
%                 end
%                 rwArray{size(rwArray, 1), 2} = [randi([1, 10]), bin2dec(strcat(dec2bin(round(7.5 + 0.5*rn)), dec2bin(randi([0, 15]))))];
%                 rwArray{size(rwArray, 1), 2} = [1, bin2dec(strcat(dec2bin(round(7.5 + 0.5*rn)), dec2bin(randi([0, 15]))))];
            else
                rwArray{size(rwArray, 1) + 1, 1} = 'r';
                rwArray{size(rwArray, 1), 2} = randi([1, numAddresses]);
            end
        end
        [writeStallTime, readStallTime, totalTime, readE, writeE] = MMU_rw_parallel_wEnergy(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult, writeBufferSize, maxWaitBufferSize, bitReadTime, switchDetectTime, writeVolt, I_LRS, I_HRS, readEnergyPerBit, bufferAccessEnergy, Clk_M_Value, 0, debug);
        times_parallel(i, j) = totalTime;
        % output = MMU_rw_parallel_new(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMultBase, 1, bitReadTime, switchDetectTime, 0, 1);
        times_parallel_base(i, j) = numAddresses * 10;
%         times_parallel(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult, writeBufferSize, 0);
%         times_serial(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult*8, writeBufferSize, 0);
%         times_parallel_base(i, j) = MMU_rw_parallel(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase, 1, 1);
%         times_serial_base(i, j) = MMU_rw_serial(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMultBase*8, 1, 1);
    end
    fprintf('%0.1f%% Completed\n', i/length(rwSize)*100);
end

figure;
errorbar(rwSize, 100.*(linspace(1, 1, size(rwSize, 2))' - mean(times_parallel, 2)./mean(times_parallel_base, 2)), 100.*(linspace(1, 1, size(rwSize, 2))' - mean(times_parallel, 2)./mean(times_parallel_base, 2)) - 100.*(linspace(1, 1, size(rwSize, 2))' - max(times_parallel, [], 2)./min(times_parallel_base, [], 2)), 100.*(linspace(1, 1, size(rwSize, 2))' - min(times_parallel, [], 2)./max(times_parallel_base, [], 2)) - 100.*(linspace(1, 1, size(rwSize, 2))' - mean(times_parallel, 2)./mean(times_parallel_base, 2)));
hold on;
plot(rwSize, 100.*(linspace(1, 1, size(rwSize, 2))' - mean(times_parallel, 2)./mean(times_parallel_base, 2)));
hold off;
xlabel('Read/Write Cluster Size');
ylabel('Performance Gain (%)');


% figure;
% errorbar(rwSize, mean(times_parallel, 2), mean(times_parallel, 2) - min(times_parallel, [], 2), max(times_parallel, [], 2) - mean(times_parallel, 2), '-o');
% xlabel('Read/Write Cluster Size');
% ylabel('Write time (Clk_Ms)');
% 
% hold on;
% errorbar(rwSize, mean(times_parallel_base, 2), mean(times_parallel_base, 2) - min(times_parallel_base, [], 2), max(times_parallel_base, [], 2) - mean(times_parallel_base, 2), '-o');
% hold off;
% legend({'Parallel Write Proposed (8 bits)', 'Parallel Write Reference (8 bits)'})
% 
% figure;
% errorbar(rwSize, mean(times_serial, 2), mean(times_serial, 2) - min(times_serial, [], 2), max(times_serial, [], 2) - mean(times_serial, 2), '-o');
% xlabel('Read/Write Cluster Size');
% ylabel('Write time (Clk_Ms)');
% 
% hold on;
% errorbar(rwSize, mean(times_serial_base, 2), mean(times_serial_base, 2) - min(times_serial_base, [], 2), max(times_serial_base, [], 2) - mean(times_serial_base, 2), '-o');
% hold off;
% legend({'Serial Write Proposed', 'Serial Write Reference'})

