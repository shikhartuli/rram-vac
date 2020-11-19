function [outputTime, readE, writeE, writeTimes_inst, totalEnergy_inst] = MMU_rw_parallel_wEnergy_wInst(rwArray, mu, sigma, numBits, leftDiscard, rightDiscard, offsetMult, writeBufferSize, bitReadTime, switchDetectTime, writeVolt, I_LRS, I_HRS, readEnergyPerBit, Clk_M_Value, base, debug)

% leftDsicard belongs from 0 to 1. If writeTime > leftDiscard * Clk_M and 
% writeTime < Clk_M, then these values truncated to Clk_M
% rightDiscard is in sigma's, i.e. rightDiscard = 3 means 3*sigma right of
% mu*Clk_M

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

Clk_M = 1;

% mu = 8; % Mean of the switching time distribution
% sigma = 3;  % Sigma of the switching time distribution
% leftDiscard = 0;
% rightDiscard = 3;
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
% Clk_M_Value = 5e-9; % in seconds
% debug = 1;

defaultWordInMemory = '0';
for n = 1:(numBits-1)
    defaultWordInMemory = strcat(defaultWordInMemory, '0');
end

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
% rwArray = {'w', [1, 1]; 'w', [2, 2]; 'r', 2; 'w', [3, 3]; 'w', [4, 4]; 'w', [2, 7]; 'w', [5, 255]; 'w', [6, 254]; 'w', [5, 255]; 'r', 1; 'w', [6, 6]};
% rwArray = {'w', [1, 10011001]; 'w', [2, 10001000]; 'w', [3, 10000000]; 'w', [4, 10101010]; 'r', 4}; 

P_cycles = size(rwArray, 1);

writeTimes_inst = zeros(P_cycles, 1);
totalEnergy_inst = zeros(P_cycles, 1);

waitBuffer = {};
% global writeBuffer;
writeBuffer = {};

% global Memory;
Memory = {};
readBuffer = {};
readMemLock = 0;
readEnergy = 0;
writeEnergy = 0;

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
    if i <= P_cycles && rwArray{i, 1} == 'w'
        if size(waitBuffer, 1) == 0
            waitBuffer(size(waitBuffer, 1) + 1, 2) = rwArray(i, 2:end);
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
        if size(Memory, 1) > 1 && time / timeNorm < Memory{end, 1}
            Memory_end = Memory{end - 1, 2};
        elseif size(Memory, 1) > 0
            Memory_end = Memory{end, 2};
        else
            Memory_end = [];
        end
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
        elseif size(Memory_end, 1) > 0 && ismember(rwArray{i, 2}, Memory_end(:, 1))
            if finishTime > Clk_P
                time = time + Clk_P;
                finishTime = finishTime - Clk_P;
                i = i - 1;
                continue;
            elseif finishTime > 0 && (Clk_P - finishTime) < bitReadTime
                i = i - 1;
                readMemLock = 1;
            elseif finishTime > 0
                for j = 1:size(Memory_end, 1)
                    if rwArray{i, 2} == Memory_end(j, 1)
                        readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  Memory_end(j, 2)];
                        readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
                        readBuffer{size(readBuffer, 1), 3} = 'from Memory';
                        readEnergy = readEnergy + numBits * readEnergyPerBit;
                        % totalEnergy_inst(i) = numBits * readEnergyPerBit;
                    end
                end
            else
                for j = 1:size(Memory_end, 1)
                    if rwArray{i, 2} == Memory_end(j, 1)
                        readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  Memory_end(j, 2)];
                        readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
                        readBuffer{size(readBuffer, 1), 3} = 'from Memory';
                        readEnergy = readEnergy + numBits * readEnergyPerBit;
                        % totalEnergy_inst(i) = numBits * readEnergyPerBit;
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
        end
    end
            
    while 1
%         if size(waitBuffer{end, 2}, 1) < writeBufferSize && finishTime >= Clk_P
%             if mod(time, Clk_P) == 0
%                 time = (time/Clk_P + 1)*Clk_P;
%                 if finishTime > 0
%                     finishTime = finishTime - Clk_P;
%                 end
%             else
%                 if finishTime > 0
%                     finishTime = finishTime - (Clk_P - time + floor(time/Clk_P)*Clk_P);
%                 end
%                 time = ceil(time/Clk_P)*Clk_P;
%             end
%             break;
%         end
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
                if size(Memory, 1) == 0
                    Memory{size(Memory, 1) + 1, 2} = writeBuffer{end, 2};
                else
                    Memory{size(Memory, 1) + 1, 2} = [Memory{size(Memory, 1), 2}; writeBuffer{end, 2}];
                end
                
                % Maintaining proper write operation for the change of data
                % at the same address
                for j = 1:size(Memory{end, 2}, 1)
                    for k = j:size(Memory{end, 2}, 1)
                        if k > j && Memory{end, 2}(j, 1) == Memory{end, 2}(k, 1)
                            Memory{end, 2}(j, 2) = Memory{end, 2}(k, 2);
                            Memory{end, 2}(k, :) = [];
                            break;
                        end
                    end
                end
                    
                Memory(size(Memory, 1), 1) = {(time + finishTime) / timeNorm};
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

if size(readBuffer, 1) > 0
    timeReadBuffer = readBuffer{end, 1};
else
    timeReadBuffer = 0;
end
if size(Memory, 1) > 0
    timeMemory = Memory{end, 1};
else
    timeMemory = 0;
end

outputTime = max(timeReadBuffer, timeMemory);
readE = readEnergy;
writeE = writeEnergy;

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
                [tf, idx] = ismember(wordAddress, Memory{end, 2}(:, 1));
                if tf
                    wordInMemory = dec2bin(Memory{end, 2}(idx, 2), numBits);
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
                        totalEnergy_inst(i) = totalEnergy_inst(i) + Clk_M_Value * writeVolt * I_HRS * Clk_P;
                    elseif bitToWrite == '1'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * I_LRS * Clk_P;
                        totalEnergy_inst(i) = totalEnergy_inst(i) + Clk_M_Value * writeVolt * I_LRS * Clk_P;
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
                        totalEnergy_inst(i) = totalEnergy_inst(i) + Clk_M_Value * writeVolt * (I_HRS * (writeTime_test - switchDetectTime) + I_LRS * (Clk_P - (writeTime_test - switchDetectTime)));
                    elseif wordInMemory(j) == '1' && bitToWrite == '0'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * (I_LRS * (writeTime_test - switchDetectTime) + I_HRS * (Clk_P - (writeTime_test - switchDetectTime)));
                        totalEnergy_inst(i) = totalEnergy_inst(i) + Clk_M_Value * writeVolt * (I_LRS * (writeTime_test - switchDetectTime) + I_HRS * (Clk_P - (writeTime_test - switchDetectTime)));
                    end
                end
            end
        end
        finishTime = size(writeBuffer{end, 2}, 1) * Clk_P; % (mu + rightDiscard * sigma) * Clk_M;
        writeTimes_inst(i) = finishTime;
%         for inst = i:(i + floor(finishTime/Clk_P))
%             writeTimes_inst(inst) = finishTime/size(writeBuffer{end, 2}, 1);
%         end
    else
        for w = 1:size(writeBuffer{end, 2}, 1)
            wordToWrite = dec2bin(writeBuffer{end, 2}(w, 2), numBits);
            wordAddress = writeBuffer{end, 2}(w, 1);
            if size(Memory, 1) > 0
                [tf, idx] = ismember(wordAddress, Memory{end, 2}(:, 1));
                if tf
                    wordInMemory = dec2bin(Memory{end, 2}(idx, 2), numBits);
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
                        totalEnergy_inst(i) = totalEnergy_inst(i) + Clk_M_Value * writeVolt * I_HRS * writeTime_test;
                    elseif bitToWrite == '1'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * I_LRS * writeTime_test;
                        totalEnergy_inst(i) = totalEnergy_inst(i) + Clk_M_Value * writeVolt * I_LRS * writeTime_test;
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
                        totalEnergy_inst(i) = totalEnergy_inst(i) + Clk_M_Value * writeVolt * (I_HRS * (writeTime_test - switchDetectTime) + I_LRS * switchDetectTime);
                    elseif wordInMemory(j) == '1' && bitToWrite == '0'
                        writeEnergy = writeEnergy + Clk_M_Value * writeVolt * (I_LRS * (writeTime_test - switchDetectTime) + I_HRS * switchDetectTime);
                        totalEnergy_inst(i) = totalEnergy_inst(i) + Clk_M_Value * writeVolt * (I_LRS * (writeTime_test - switchDetectTime) + I_HRS * switchDetectTime);
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
        writeTimes_inst(i) = finishTime;
%         for inst = i:(i + floor(finishTime/Clk_P))
%             writeTimes_inst(inst) = finishTime/size(writeBuffer{end, 2}, 1);
%         end
    end

    end

end

