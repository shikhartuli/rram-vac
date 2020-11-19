% function output = MMU_rw(rwArray, mu, sigma, leftDiscard, rightDiscard, offsetMult, writeBufferSize)

% This code assumes that the write/read operations are one with only
% one-bit data

% leftDsicard belongs from 0 to 1. If writeTime > leftDiscard * Clk_M and 
% writeTime < Clk_M, then these values truncated to Clk_M
% rightDiscard is in sigma's, i.e. rightDiscard = 3 means 3*sigma right of
% mu*Clk_M

mu = 8; % Mean of the switching time distribution
sigma = 3;  % Sigma of the switching time distribution
leftDiscard = 0;
rightDiscard = 3;
writeBufferSize = 2;
offsetMult = 1.1;

timeDisplay = 'Clk_P';
timeNorm = 1;

if timeDisplay == 'Clk_P'
    timeNorm = offsetMult * mu;
else
    timeNorm = 1;
end

Clk_M = 1;
Clk_P = offsetMult * mu * Clk_M; % Clk_P = mu * Clk_M for no accumulation 
% in steady state

rwArray = {'w', [1, 1]; 'w', [2, 0]; 'w', [3, 1]; 'w', [4, 0]; 'r', 4};

P_cycles = size(rwArray, 1);

waitBuffer = {};
global writeBuffer;
writeBuffer = {};
Memory = {};
readBuffer = {};
fillReadBuffer = 0;

overflow = 0; % Overflow Memory cycles due to asynchronous nature
global time;
time = 0;
% writeTime = [];
finishTime = 0;
global finishTimeArray;
finishTimeArray = [];

i = 0;

while i <= P_cycles+2
    i = i + 1;
    global time;
    if i <= P_cycles && rwArray{i, 1} == 'w'
        if size(waitBuffer, 1) == 0
            waitBuffer(size(waitBuffer, 1) + 1, 2) = rwArray(i, 2:end);
        else
            % CAM write operations maintins only one data for a particular
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
        % But, the memory and the writeBuffer are filled ant emptied at 
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
        if size(waitBuffer{end, 2}, 1) > 0 && ismember(rwArray{i, 2}, waitBuffer{end, 2}(:, 1))
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
                i = i - 1;
                continue;
            else
                for j = 1:size(Memory_end, 1)
                    if rwArray{i, 2} == Memory_end(j, 1)
                        readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  Memory_end(j, 2)];
                        readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
                        readBuffer{size(readBuffer, 1), 3} = 'from Memory';
                    end
                end
            end
        else
            readBuffer{size(readBuffer, 1) + 1, 2} = [rwArray{i, 2},  NaN];
            readBuffer(size(readBuffer, 1), 1) = {(time + Clk_P) / timeNorm};
            readBuffer{size(readBuffer, 1), 3} = 'Address data not found';
        end
    end
            
    while 1
        if size(waitBuffer{end, 2}, 1) < writeBufferSize && finishTime >= Clk_P
            if mod(time, Clk_P) == 0
                time = (time/Clk_P + 1)*Clk_P;
                if finishTime > 0
                    finishTime = finishTime - Clk_P;
                end
            else
                if finishTime > 0
                    finishTime = finishTime - (Clk_P - time + floor(time/Clk_P)*Clk_P);
                end
                time = ceil(time/Clk_P)*Clk_P;
            end
            break;
        end
%         if finishTime >= Clk_P
%             if mod(time, Clk_P) == 0
%                 time = (time/Clk_P + 1)*Clk_P;
%                 finishTime = finishTime - Clk_P;
%             else
%                 finishTime = finishTime - (Clk_P - time + floor(time/Clk_P)*Clk_P);
%                 time = ceil(time/Clk_P)*Clk_P;
%             end
%             break;
%         end
        if finishTime >= 0
            % time = time + finishTime; % Important assumption that shifing 
            % operation is instantaneous
            
            if size(writeBuffer, 1) > 0 && size(writeBuffer{end, 2}, 1) > 0
                if size(Memory, 1) == 0
                    Memory{size(Memory, 1) + 1, 2} = writeBuffer{end, 2};
                else
                    Memory{size(Memory, 1) + 1, 2} = [Memory{size(Memory, 1), 2}; writeBuffer{end, 2}];
                end
                Memory(size(Memory, 1), 1) = {(time + finishTime) / timeNorm};
                writeBuffer{size(writeBuffer, 1) + 1, 2} = [];
                writeBuffer(size(writeBuffer, 1), 1) = {(time + finishTime) / timeNorm};
            end
            
            if size(waitBuffer{end, 2}, 1) >= writeBufferSize
                writeBuffer{size(writeBuffer, 1) + 1, 2} = waitBuffer{end, 2}(1:writeBufferSize, :);
                writeBuffer(size(writeBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                waitBuffer{size(waitBuffer, 1) + 1, 2} = waitBuffer{end, 2}(writeBufferSize + 1:end, :);
                waitBuffer(size(waitBuffer, 1), 1) = {(time + finishTime) / timeNorm};
                
                finishTime = getFinishTime(writeBufferSize, Clk_M, Clk_P, leftDiscard, rightDiscard, mu, sigma) ...
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

% end

function finishTime = getFinishTime(writeBufferSize, Clk_M, Clk_P, leftDiscard, rightDiscard, mu, sigma)

global time;
global writeBuffer;
global finishTimeArray;
finishTime = 0;
writeTime = [];

for i = 1:writeBufferSize
    writeTime_test = Clk_M * round(mu + sigma*randn);
    while writeTime_test > (mu + rightDiscard * sigma) * Clk_M || writeTime_test < leftDiscard * Clk_M
        writeTime_test = Clk_M * round(mu + sigma*randn);
    end
    writeTime = [writeTime; writeTime_test];   % > Right side discarded
    if writeTime(end) > leftDiscard * Clk_M && writeTime(end) < Clk_M
        writeTime(end) = Clk_M;
    end
    finishTime = writeTime(end) + finishTime;
    finishTimeArray = [finishTimeArray; writeTime(end)];
end

end
