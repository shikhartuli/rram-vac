function output = MMU_rw_parallel_new(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult, writeBufferSize, bitReadTime, switchDetectTime, maxWaitBufferSize, base)

warning('off', 'all');
% leftDsicard belongs from 0 to 1. If writeTime > leftDiscard * Clk_M and 
% writeTime < Clk_M, then these values truncated to Clk_M
% rightDiscard is in sigma's, i.e. rightDiscard = 3 means 3*sigma right of
% mu*Clk_M

Clk_M = 1;
global numBits;

% mu = 5; % Mean of the switching time distribution
% sigma = 2;  % Sigma of the switching time distribution
% leftDiscard = 0;
% rightDiscard = 2;
% writeBufferSize = 6;
% offsetMult = 1;
% base = 0;
% numBits = 8;
% bitReadTime = 1* Clk_M;
% switchDetectTime = 0.2 * Clk_M;
% maxWaitBufferSize = 25;
debug = 0;

timeDisplay = 'Clk_M';
timeNorm = 1;

if timeDisplay == 'Clk_P'
    timeNorm = offsetMult * mu;
else
    timeNorm = 1;
end

Clk_P = offsetMult * mu * Clk_M; % Clk_P = mu * Clk_M for no accumulation 
% in steady state

% rwArray = {'w', [1, 0]; 'w', [2, 0]; 'w', [3, 0]; 'w', [1, 255]; 'r', 1; 'r', 3};
% rwArray = {'w', [1, 1]; 'w', [2, 2]; 'w', [3, 3]; 'w', [4, 4]; 'w', [2, 7]; 'w', [5, 255]; 'w', [6, 254]; 'w', [5, 255]; 'r', 1; 'w', [6, 6]};
% rwArray = {'w', [1, 10011001]; 'w', [2, 10001000]; 'w', [3, 10000000]; 'w', [4, 10101010]; 'r', 4}; 

P_cycles = size(rwArray, 1);

waitBuffer = {};
global writeBuffer;
writeBuffer = {};

% global writeTime;

maxAddress = 0;
for i = 1:1:P_cycles
    if rwArray{i, 2}(1, 1) > maxAddress
        maxAddress = rwArray{i, 2}(1, 1);
    end
end

Memory = NaN(maxAddress, 1);
% global Memory;
% Memory = {};
readBuffer = {};
readMemLock = 0;
readEnergy = 0;
writeEnergy = 0;
writeStallCycles = 0;
readStallCycles = 0;
flushWrite = 0;

overflow = 0; % Overflow Memory cycles due to asynchronous nature
% global time;
time = 0;
% writeTime = [];
finishTime = 0;
% global finishTimeArray;
% finishTimeArray = [];

i = 0;

while i <= P_cycles+(4*writeBufferSize*8)
    if debug == 0 && size(writeBuffer, 1) > 10
        writeBuffer = writeBuffer(5:end, :);
    end
    if debug == 0 && size(waitBuffer, 1) > 10
        waitBuffer = waitBuffer(5:end, :);
    end
    if debug == 0 && size(readBuffer, 1) > 10
        readBuffer = readBuffer(5:end, :);
    end
    i = i + 1;
    global time;
    if i <= P_cycles && rwArray{i, 1} == 'w'
        if size(waitBuffer, 1) == 0
            waitBuffer(size(waitBuffer, 1) + 1, 2) = rwArray(i, 2:end);
        elseif maxWaitBufferSize > 0 && finishTime == 0 && size(waitBuffer{end, 2}, 1) == maxWaitBufferSize
            writeStallCycles = writeStallCycles + 1;
            i = i - 1;
        elseif maxWaitBufferSize > 0 && finishTime > 0 && size(waitBuffer{end, 2}, 1) == (maxWaitBufferSize - writeBufferSize)
            writeStallCycles = writeStallCycles + 1;
            i = i - 1;
        else
            % CAM write operations maintains only one data for a particular
            % address in the waitBuffer
            if ismember(rwArray{i, 2}(1, 1), waitBuffer{end, 2}(:, 1))
                [~, idx] = ismember(rwArray{i, 2}(1, 1), waitBuffer{end, 2}(:, 1));
                waitBuffer{size(waitBuffer, 1) + 1, 2} = waitBuffer{end, 2};
                waitBuffer{size(waitBuffer, 1), 2}(idx, 2) = rwArray{i, 2}(1, 2);
            else
                waitBuffer{size(waitBuffer, 1) + 1, 2} = [waitBuffer{end, 2}; rwArray{i, 2:end}];
            end
            if i == P_cycles
                if mod(time, Clk_P) == 0
                    flushWrite = (time/Clk_P + 1)*Clk_P;
                else
                    flushWrite = ceil(time/Clk_P)*Clk_P;
                end
            end
        end
        waitBuffer(size(waitBuffer, 1), 1) = {time / timeNorm};
    elseif i <= P_cycles && rwArray{i, 1} == 'r'
        % At this point of time is a multiple of Clk_P
        % waitBuffer is always filled at multiples of Clk_P (not emptied)
        % But, the memory and the writeBuffer are filled and emptied at 
        % arbitrary times. So we need to check that the current time 
        % is < the Memory{end, 1} or writeBuffer{end, 1}
        if size(writeBuffer, 1) > 1 && time / timeNorm < writeBuffer{end, 1} 
            writeBuffer_end = writeBuffer{end - 1, 2};
        elseif size(writeBuffer, 1) > 0
            writeBuffer_end = writeBuffer{end, 2};
        else 
            writeBuffer_end = [];
        end
%         if size(Memory, 1) > 1 && time / timeNorm < Memory{end, 1}
%             Memory_end = Memory{end - 1, 2};
%         elseif size(Memory, 1) > 0
%             Memory_end = Memory{end, 2};
%         else
%             Memory_end = [];
%         end
        if size(waitBuffer, 1) > 0 && size(waitBuffer{end, 2}, 1) > 0 && ismember(rwArray{i, 2}, waitBuffer{end, 2}(:, 1))
            for j = 1:size(waitBuffer{end, 2}, 1)
                if rwArray{i, 2} == waitBuffer{end, 2}(j, 1)
                    readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2}, waitBuffer{end, 2}(j, 2)];
                    readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
                    readBuffer{size(readBuffer, 1), 3} = 'from waitBuffer';
                end
            end
        elseif size(writeBuffer_end, 1) > 0 && ismember(rwArray{i, 2}, writeBuffer_end(:, 1))
            for j = 1:size(writeBuffer_end, 1)
                if rwArray{i, 2} == writeBuffer_end(j, 1)
                    readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2}, writeBuffer_end(j, 2)];
                    readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
                    readBuffer{size(readBuffer, 1), 3} = 'from writeBuffer';
                end
            end    
        elseif ~isnan(Memory(rwArray{i, 2})) % size(Memory_end, 1) > 0 && ismember(rwArray{i, 2}, Memory_end(:, 1))
            if finishTime > Clk_P
                time = time + Clk_P;
                finishTime = finishTime - Clk_P;
                readStallCycles = readStallCycles + 1;
                i = i - 1;
                continue;
            elseif finishTime > 0 && (Clk_P - finishTime) < bitReadTime
                i = i - 1;
                readStallCycles = readStallCycles + 1;
                readMemLock = 1;
            elseif finishTime > 0
                readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  Memory(rwArray{i, 2})];
                readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
                readBuffer{size(readBuffer, 1), 3} = 'from Memory';
%                 for j = 1:size(Memory_end, 1)
%                     if rwArray{i, 2} == Memory_end(j, 1)
%                         readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  Memory_end(j, 2)];
%                         readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
%                         readBuffer{size(readBuffer, 1), 3} = 'from Memory';
%                     end
%                 end
            else
                readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  Memory(rwArray{i, 2})];
                readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
                readBuffer{size(readBuffer, 1), 3} = 'from Memory';
%                 for j = 1:size(Memory_end, 1)
%                     if rwArray{i, 2} == Memory_end(j, 1)
%                         readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  Memory_end(j, 2)];
%                         readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
%                         readBuffer{size(readBuffer, 1), 3} = 'from Memory';
%                     end
%                 end
                if mod(time, Clk_P) == 0
                    time = (time/Clk_P + 1)*Clk_P;
                else
                    time = ceil(time/Clk_P)*Clk_P;
                end
                continue;
            end
        else
            readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  NaN];
            readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
            readBuffer{size(readBuffer, 1), 3} = 'Address data not found';
        end
    end
            
    while 1
        if finishTime >= Clk_P
            if mod(time, Clk_P) == 0
                time = (time/Clk_P + 1)*Clk_P;
                finishTime = finishTime - Clk_P;
            else
                finishTime = finishTime - (Clk_P - time + floor(time/Clk_P)*Clk_P);
                time = ceil(time/Clk_P)*Clk_P;
            end
            break;
        end
        if finishTime >= 0
            % time = time + finishTime; % Important assumption that shifing 
            % operation is instantaneous
            
            if size(writeBuffer, 1) > 0 && size(writeBuffer{end, 2}, 1) > 0
%                 if size(Memory, 1) == 0
%                     Memory{size(Memory, 1) + 1, 2} = writeBuffer{end, 2};
%                 else
%                     Memory{size(Memory, 1) + 1, 2} = [Memory{size(Memory, 1), 2}; writeBuffer{end, 2}];
%                 end
%                 
%                 % Maintaining proper write operation for the change of data
%                 % at the same address
%                 for j = 1:size(Memory{end, 2}, 1)
%                     for k = j:size(Memory{end, 2}, 1)
%                         if k > j && Memory{end, 2}(j, 1) == Memory{end, 2}(k, 1)
%                             Memory{end, 2}(j, 2) = Memory{end, 2}(k, 2);
%                             Memory{end, 2}(k, :) = [];
%                             break;
%                         end
%                     end
%                 end
                
                for w = 1:1:size(writeBuffer{end, 2}, 1)
                    Memory(writeBuffer{end, 2}(w, 1)) = writeBuffer{end, 2}(w, 2);
                end
                    
%                 Memory(size(Memory, 1), 1) = {(time + finishTime) / timeNorm};
                writeBuffer{size(writeBuffer, 1) + 1, 2} = [];
                writeBuffer(size(writeBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                if readMemLock == 1
                    readMemLock = 0;
                    finishTime = 0;
                    if mod(time, Clk_P) == 0
                        time = (time/Clk_P + 1)*Clk_P;
                    else
                        time = ceil(time/Clk_P)*Clk_P;
                    end
                    break;
                end
                    
            end
            
            if size(waitBuffer, 1) > 0 && size(waitBuffer{end, 2}, 1) >= writeBufferSize
                writeBuffer{size(writeBuffer, 1) + 1, 2} = waitBuffer{end, 2}(1:writeBufferSize, :);
                writeBuffer(size(writeBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                waitBuffer{size(waitBuffer, 1) + 1, 2} = waitBuffer{end, 2}(writeBufferSize + 1:end, :);
                waitBuffer(size(waitBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                finishTime = getFinishTime(writeBufferSize, Clk_M, Clk_P, leftDiscard, rightDiscard, mu, sigma, base) ...
                    + finishTime;
                
            % Clearing out remains in the waitBuffer after last Processor
            % instruction
            elseif i >= P_cycles && size(waitBuffer, 1) > 0 && size(waitBuffer{end, 2}, 1) > 0
                writeBuffer{size(writeBuffer, 1) + 1, 2} = waitBuffer{end, 2}(1:end, :);
                writeBuffer(size(writeBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                % waitBuffer{size(waitBuffer, 1) + 1, 2} = waitBuffer{end, 2}(writeBufferSize + 1:end, :);
                waitBuffer{size(waitBuffer, 1) + 1, 2} = [];
                waitBuffer(size(waitBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                finishTime = getFinishTime(writeBufferSize, Clk_M, Clk_P, leftDiscard, rightDiscard, mu, sigma, base) ...
                    + finishTime;
            else
                finishTime = 0;
                if mod(time, Clk_P) == 0
                    time = (time/Clk_P + 1)*Clk_P;
                else
                    time = ceil(time/Clk_P)*Clk_P;
                end
                break;
            end
        end
    end
end

% maxBufferSize = 0;
% for w = 1:1:size(waitBuffer, 1)
%     if size(waitBuffer{w, 2}, 1) > maxBufferSize
%         maxBufferSize = size(waitBuffer{w, 2}, 1);
%     end
% end
% output = maxBufferSize;
% output = Memory{end, 1};
% output = writeStallCycles;
% if flushWrite == 0
%     output = [writeStallCycles, Memory{end, 1}];
% else
%     output = [writeStallCycles, flushWrite];
% end
output = [writeStallCycles, flushWrite];
end

function finishTime = getFinishTime(writeBufferSize, Clk_M, Clk_P, leftDiscard, rightDiscard, mu, sigma, base)

global time;
global writeBuffer;
global Memory;
global numBits;
finishTime = 0;
global writeTime;
writeTime = zeros(writeBufferSize, numBits);

if base == 1
    finishTime = size(writeBuffer{end, 2}, 1) * Clk_P; % (mu + rightDiscard * sigma) * Clk_M;
    
else
    for i = 1:size(writeBuffer{end, 2}, 1)
        wordToWrite = dec2bin(writeBuffer{end, 2}(i, 2), numBits);
        wordAddress = writeBuffer{end, 2}(i, 1);
        for j = 1:numBits
            bitToWrite = wordToWrite(j);
            if size(Memory, 1) > 0
%                 [tf, idx] = ismember(wordAddress, Memory{end, 2}(:, 1));
                tf = ~isnan(Memory(wordAddress));
                if tf
                    wordInMemory = dec2bin(Memory(wordAddress), numBits);
                end
            else
                tf = false(1);
            end
            if tf && bitToWrite == wordInMemory(j)
                writeTime_test = Clk_M;
            else
                writeTime_test = Clk_M * round(mu + sigma*randn);
                while writeTime_test > (mu + rightDiscard * sigma) * Clk_M || writeTime_test < leftDiscard * Clk_M
                    writeTime_test = Clk_M * round(mu + sigma*randn);
                end
                if writeTime_test > leftDiscard * Clk_M && writeTime_test < Clk_M
                    writeTime_test = Clk_M;
                end
            end
            if i > 1
                writeTime(i, j) = writeTime(i-1, j) + writeTime_test;
            else
                writeTime(i, j) = writeTime_test;
            end
        end
    end
    finishTime = max(writeTime(size(writeBuffer{end, 2}, 1), :));
end

end
