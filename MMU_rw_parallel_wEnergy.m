function [writeStallTime, readStallTime, totalTime, readE, writeE] = MMU_rw_parallel_wEnergy(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult, writeBufferSize, waitBufferSize, bitReadTime, switchDetectTime, writeVolt, I_LRS, I_HRS, readEnergyPerBit, bufferAccessEnergy, Clk_M_Value, base, debug)

% leftDiscard belongs from 0 to 1. If writeTime > leftDiscard * Clk_M and 
% writeTime < Clk_M, then these values truncated to Clk_M
% rightDiscard is in sigma's, i.e. rightDiscard = 3 means 3*sigma right of
% mu*Clk_M

Clk_M = 1;

% Only un-comment the lines below while debugging
% global numBits;
% global switchDetectTime;
% global readEnergy;
% global writeEnergy;
% global I_HRS;
% global I_LRS;
% global writeVolt;
% global Clk_M_Value;
% global writeBufferSize;
% global base;

% mu = 8; % Mean of the switching time distribution
% sigma = 3;  % Sigma of the switching time distribution
% leftDiscard = 0;
% rightDiscard = 3;
% waitBufferSize = 3;
% writeBufferSize = 1;
% offsetMult = 1.1;
% base = 0;
% numBits = 32;
% bitReadTime = 1* Clk_M;
% switchDetectTime = 0.2 * Clk_M;
% writeVolt = 1; % in Voltts
% I_LRS = 100e-6; % in Amps
% I_HRS = 15e-6; % in Amps
% readEnergyPerBit = 1e-12; % in J
% bufferAccessEnergy = 200e-15;
% Clk_M_Value = 5e-9; % in seconds
% debug = 1;

defaultWordInMemory = '0';
for n = 1:(numBits-1)
    defaultWordInMemory = strcat(defaultWordInMemory, '0');
end

timeDisplay = 'Clk_M';  % Alters display output in terms of Clk_P or
timeNorm = 1;           % Clk_M. Keep this as Clk_M, unless debugging

if timeDisplay == 'Clk_P'
    timeNorm = offsetMult * mu;
else
    timeNorm = 1;
end

Clk_P = offsetMult * mu * Clk_M; % Clk_P = mu * Clk_M for no accumulation 
                                 % in steady state

% Only un-comment the lines below while debugging
% rwArray = {'w', [1, 0]; 'w', [2, 0]; 'w', [3, 0]; 'w', [1, 255]; 'r', 1; 'r', 3};
% rwArray = {'w', [1, 1]; 'w', [2, 2]; 'r', 2; 'w', [3, 3]; 'w', [4, 4]; 'w', [2, 7]; 'w', [5, 255]; 'w', [6, 254]; 'w', [5, 255]; 'r', 1; 'w', [6, 6]};
% rwArray = {'w', [1, 10011001]; 'w', [2, 10001000]; 'w', [3, 10000000]; 'w', [4, 10101010]; 'r', 4}; 

P_cycles = size(rwArray, 1);

% global waitBuffer;
waitBuffer = {};
% global writeBuffer;
writeBuffer = {};

maxAddress = 0;
for i = 1:1:P_cycles
    if rwArray{i, 2}(1, 1) > maxAddress
        maxAddress = rwArray{i, 2}(1, 1);
    end
end

% Initializing Memory with the size equal to the maximum address in the
% instructions
Memory = NaN(maxAddress, 1);
% global Memory;
% Memory = {};
readBuffer = {};

% readMemLock is 1 when the Memory can't be accessed for reads due to
% "Lock", while being written
readMemLock = 0;
readEnergy = 0;
writeEnergy = 0;

overflow = 0; % Overflow Memory cycles due to asynchronous nature
% global time;
time = 0;
finishTime = 0;
flushWrite = 0;
writeStallCycles = 0;
readStallCycles = 0;

i = 0;

while i <= P_cycles+(4*writeBufferSize*8)
    
    % When "debug" is 0, only the last few entries of the writeBuffer, 
    % readBuffer and the Memory are maintained
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
    if i <= P_cycles && rwArray{i, 1} == 'w'
        if base == 0
            % Each time the buffer is accessed, an extra energy is spent
            writeEnergy = writeEnergy + bufferAccessEnergy;
        end
        if size(waitBuffer, 1) == 0
            waitBuffer(size(waitBuffer, 1) + 1, 2) = rwArray(i, 2:end);
        % If the writeBuffer is "Locked", i.e. finishTime > 0, then the 
        % waitBuffer can only increase till (waitBufferSize -
        % writeBufferSize), otherwise, it can increase till waitBufferSize
        elseif waitBufferSize > 0 && finishTime == 0 && size(waitBuffer{end, 2}, 1) == waitBufferSize
            writeStallCycles = writeStallCycles + 1;
            i = i - 1;
        elseif waitBufferSize > 0 && finishTime > 0 && size(waitBuffer{end, 2}, 1) == (waitBufferSize - writeBufferSize)
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
            % flushWrite is the time  when all processor operations are
            % complete and only flushing operations from the buffer to the
            % Memory are left
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
            if base == 0
                % Again, when the buffer is read, we need to add the
                % bufferAccessEnergy
                readEnergy = readEnergy + bufferAccessEnergy;
            end
            if i == P_cycles
                if mod(time, Clk_P) == 0
                    flushWrite = (time/Clk_P + 1)*Clk_P;
                else
                    flushWrite = ceil(time/Clk_P)*Clk_P;
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
            if base == 0
                readEnergy = readEnergy + bufferAccessEnergy;
            end
            if i == P_cycles
                if mod(time, Clk_P) == 0
                    flushWrite = (time/Clk_P + 1)*Clk_P;
                else
                    flushWrite = ceil(time/Clk_P)*Clk_P;
                end
            end
        elseif ~isnan(Memory(rwArray{i, 2})) % size(Memory_end, 1) > 0 && ismember(rwArray{i, 2}, Memory_end(:, 1))
            if finishTime > Clk_P
                time = time + Clk_P;
                finishTime = finishTime - Clk_P;
                i = i - 1;
                % readStallCycles represents the number of cycles the
                % processot has to stall, while the Memory is "Locked" for
                % writing operations
                readStallCycles = readStallCycles + 1;
                continue;
            elseif finishTime > 0 && (Clk_P - finishTime) < bitReadTime
                i = i - 1;
                readStallCycles = readStallCycles + 1;
                readMemLock = 1;
            elseif finishTime > 0
                readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  Memory(rwArray{i, 2})];
                readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
                readBuffer{size(readBuffer, 1), 3} = 'from Memory';
                
                readEnergy = readEnergy + numBits * readEnergyPerBit;
%                 for j = 1:size(Memory_end, 1)
%                     if rwArray{i, 2} == Memory_end(j, 1)
%                         readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  Memory_end(j, 2)];
%                         readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
%                         readBuffer{size(readBuffer, 1), 3} = 'from Memory';
%                         readEnergy = readEnergy + numBits * readEnergyPerBit;
%                     end
%                 end
                if base == 0
                    readEnergy = readEnergy + bufferAccessEnergy;
                end
                if i == P_cycles
                    if mod(time, Clk_P) == 0
                        flushWrite = (time/Clk_P + 1)*Clk_P;
                    else
                        flushWrite = ceil(time/Clk_P)*Clk_P;
                    end
                end
            else
                readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  Memory(rwArray{i, 2})];
                readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
                readBuffer{size(readBuffer, 1), 3} = 'from Memory';
                
                readEnergy = readEnergy + numBits * readEnergyPerBit;
%                 for j = 1:size(Memory_end, 1)
%                     if rwArray{i, 2} == Memory_end(j, 1)
%                         readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  Memory_end(j, 2)];
%                         readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
%                         readBuffer{size(readBuffer, 1), 3} = 'from Memory';
%                         readEnergy = readEnergy + numBits * readEnergyPerBit;
%                     end
%                 end
                if base == 0
                    readEnergy = readEnergy + bufferAccessEnergy;
                end
                if i == P_cycles
                    if mod(time, Clk_P) == 0
                        flushWrite = (time/Clk_P + 1)*Clk_P;
                    else
                        flushWrite = ceil(time/Clk_P)*Clk_P;
                    end
                end
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
            if base == 0
                readEnergy = readEnergy + bufferAccessEnergy;
            end
            if i == P_cycles
                if mod(time, Clk_P) == 0
                    flushWrite = (time/Clk_P + 1)*Clk_P;
                else
                    flushWrite = ceil(time/Clk_P)*Clk_P;
                end
            end
        end
    end
    
    % This while loop takes into account all the transfers between the
    % waitBuffer and the writeBuffer, between the writeBuffer and the 
    % Memory in between the cycles
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
            % Important assumption that shifing 
            % operation is instantaneous since in real implementation, only
            % a "Lock" bit is switched
            
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
%                     
%                 Memory(size(Memory, 1), 1) = {(time + finishTime) / timeNorm};

                for w = 1:1:size(writeBuffer{end, 2}, 1)
                    Memory(writeBuffer{end, 2}(w, 1)) = writeBuffer{end, 2}(w, 2);
                end

                writeBuffer{size(writeBuffer, 1) + 1, 2} = [];
                writeBuffer(size(writeBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                % Removes the readMemLock once the Memory write operations
                % are finished
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
            
            % Transfer from waitBuffer to the writeBuffer and determining
            % the next "finishTime" for Memory write operation
            if size(waitBuffer, 1) > 0 && size(waitBuffer{end, 2}, 1) >= writeBufferSize
                writeBuffer{size(writeBuffer, 1) + 1, 2} = waitBuffer{end, 2}(1:writeBufferSize, :);
                writeBuffer(size(writeBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                waitBuffer{size(waitBuffer, 1) + 1, 2} = waitBuffer{end, 2}(writeBufferSize + 1:end, :);
                waitBuffer(size(waitBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                finishTime = getFinishTime() ...
                    + finishTime;
                
            % Clearing out remains in the waitBuffer after last Processor
            % instruction
            elseif i >= P_cycles && size(waitBuffer, 1) > 0 && size(waitBuffer{end, 2}, 1) > 0
                writeBuffer{size(writeBuffer, 1) + 1, 2} = waitBuffer{end, 2}(1:end, :);
                writeBuffer(size(writeBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                % waitBuffer{size(waitBuffer, 1) + 1, 2} = waitBuffer{end, 2}(writeBufferSize + 1:end, :);
                waitBuffer{size(waitBuffer, 1) + 1, 2} = [];
                waitBuffer(size(waitBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                finishTime = getFinishTime() ...
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

% if size(readBuffer, 1) > 0
%     timeReadBuffer = readBuffer{end, 1};
% else
%     timeReadBuffer = 0;
% end
% if size(Memory, 1) > 0
%     timeMemory = Memory{end, 1};
% else
%     timeMemory = 0;
% end

writeStallTime = writeStallCycles;
readStallTime = readStallCycles;
totalTime = flushWrite;
readE = readEnergy;
writeE = writeEnergy;

    % This function determines the time it takes for a batch to write,
    % based on the paralle write method. When base == 1, the write is at
    % the worst case (constant) frequency, otherwise, it is according to
    % the Normal Distribution
    function finishTime = getFinishTime()

%     global writeBuffer;
%     global Memory;
%     global numBits;
%     global switchDetectTime;
%     global writeEnergy;
%     global I_HRS;
%     global I_LRS;
%     global writeVolt;
%     global Clk_M_Value;
%     global writeBufferSize;
%     global writeTime;
    finishTime = 0;
    writeTime = zeros(writeBufferSize, numBits);

    if base == 1
        for w = 1:size(writeBuffer{end, 2}, 1)
            wordToWrite = dec2bin(writeBuffer{end, 2}(w, 2), numBits);
            wordAddress = writeBuffer{end, 2}(w, 1);
            if size(Memory, 1) > 0
%                 [tf, idx] = ismember(wordAddress, Memory{end, 2}(:, 1));
                tf = ~isnan(Memory(wordAddress));
                if tf
                    wordInMemory = dec2bin(Memory(wordAddress), numBits);
                else
                    wordInMemory = defaultWordInMemory;
                end
            else
                tf = false(1);
                wordInMemory = defaultWordInMemory;
            end
            for j = 1:numBits
                bitToWrite = wordToWrite(j);
                if bitToWrite == wordInMemory(j)
                    if bitToWrite == '0'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * I_HRS * Clk_P;
                    elseif bitToWrite == '1'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * I_LRS * Clk_P;
                    end
                else
                    writeTime_test = Clk_M * round(mu + sigma*randn);
                    while writeTime_test > (mu + rightDiscard * sigma) * Clk_M || writeTime_test < leftDiscard * Clk_M
                        writeTime_test = Clk_M * round(mu + sigma*randn);
                    end
                    if writeTime_test > leftDiscard * Clk_M && writeTime_test < Clk_M
                        writeTime_test = Clk_M;
                    end
                    if wordInMemory(j) == '0' && bitToWrite == '1'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * (I_HRS * (writeTime_test - switchDetectTime) + I_LRS * (Clk_P - (writeTime_test - switchDetectTime)));
                    elseif wordInMemory(j) == '1' && bitToWrite == '0'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * (I_LRS * (writeTime_test - switchDetectTime) + I_HRS * (Clk_P - (writeTime_test - switchDetectTime)));
                    end
                end
            end
        end
        finishTime = size(writeBuffer{end, 2}, 1) * Clk_P; % (mu + rightDiscard * sigma) * Clk_M;
    else
        for w = 1:size(writeBuffer{end, 2}, 1)
            wordToWrite = dec2bin(writeBuffer{end, 2}(w, 2), numBits);
            wordAddress = writeBuffer{end, 2}(w, 1);
            if size(Memory, 1) > 0
%                 [tf, idx] = ismember(wordAddress, Memory{end, 2}(:, 1));
                tf = ~isnan(Memory(wordAddress));
                if tf
                    wordInMemory = dec2bin(Memory(wordAddress), numBits);
                else
                    wordInMemory = defaultWordInMemory;
                end
            else
                tf = false(1);
                wordInMemory = defaultWordInMemory;
            end
            for j = 1:numBits
                bitToWrite = wordToWrite(j);
                if bitToWrite == wordInMemory(j)
                    writeTime_test = switchDetectTime;
                    if bitToWrite == '0'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * I_HRS * writeTime_test;
                    elseif bitToWrite == '1'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * I_LRS * writeTime_test;
                    end
                else
                    writeTime_test = Clk_M * round(mu + sigma*randn);
                    while writeTime_test > (mu + rightDiscard * sigma) * Clk_M || writeTime_test < leftDiscard * Clk_M
                        writeTime_test = Clk_M * round(mu + sigma*randn);
                    end
                    if writeTime_test > leftDiscard * Clk_M && writeTime_test < Clk_M
                        writeTime_test = Clk_M;
                    end
                    if wordInMemory(j) == '0' && bitToWrite == '1'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * (I_HRS * (writeTime_test - switchDetectTime) + I_LRS * switchDetectTime);
                    elseif wordInMemory(j) == '1' && bitToWrite == '0'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * (I_LRS * (writeTime_test - switchDetectTime) + I_HRS * switchDetectTime);
                    end
                end
                if w > 1
                    writeTime(w, j) = writeTime(w-1, j) + writeTime_test;
                else
                    writeTime(w, j) = writeTime_test;
                end
            end
        end
        finishTime = max(writeTime(size(writeBuffer{end, 2}, 1), :));
    end

    end

end

