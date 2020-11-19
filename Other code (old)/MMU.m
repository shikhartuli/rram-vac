function output = MMU(P_cycles, mu, sigma, leftDiscard, rightDiscard, offsetMult)

% leftTruncate belongs from 0 to 1. If writeTime > leftDiscard * Clk_M and 
% writeTime < Clk_M, then these values truncated to Clk_M
% rightDiscard is in sigma's, i.e. rightDiscard = 3 means 3*sigma right of
% mu*Clk_M

% P_cycles = 20000;

% mu = 8; % Mean of the switching time distribution
% sigma = 3;  % Sigma of the switching time distribution

Clk_M = 1;
Clk_P = offsetMult * mu * Clk_M; % Clk_P = mu * Clk_M for no accumulation 
% in steady state

overflow = 0; % Overflow Memory cycles due to asynchronous nature
buffer = 0;
time = 0;
M_cycles = 0;
buffer_with_time = [];
writeTime = [];
writeToMemory= 0;

for i = 1:P_cycles
    buffer = buffer + 1;
    M_cycles = overflow; % Number of Memory cycles per Processor cycle
    buffer_with_time = [buffer_with_time; time buffer M_cycles overflow];
    while 1
        if buffer == 0
            if mod(time, Clk_P) == 0
                time = (time/Clk_P + 1)*Clk_P;
            end
            break;
        end
        if M_cycles >= Clk_P
            if mod(time, Clk_P) == 0
                time = (time/Clk_P + 1)*Clk_P;
            end
            overflow = M_cycles - Clk_P;
            break;
        end
        if overflow > 0
            buffer = buffer - 1;
            writeToMemory = writeToMemory + 1;
            time = time + overflow;
            overflow = 0;
            buffer_with_time = [buffer_with_time; time buffer M_cycles overflow];
        end
%             writeTime = [writeTime; Clk_M * round(randi([1, 2*mu-1],1,1))];
        writeTime_test = Clk_M * round(mu + sigma*randn);
        while writeTime_test > (mu + rightDiscard * sigma) * Clk_M || writeTime_test < leftDiscard * Clk_M
            writeTime_test = Clk_M * round(mu + sigma*randn);
        end
        writeTime = [writeTime; writeTime_test];   % > Right side discarded
        if writeTime(end) > leftDiscard * Clk_M && writeTime(end) < Clk_M
            writeTime(end) = Clk_M;
        end
        M_cycles = writeTime(end) + M_cycles;
        if M_cycles > Clk_P
            overflow = M_cycles - Clk_P;
            if mod(time, Clk_P) == 0
                time = (time/Clk_P + 1)*Clk_P;
            end
            break;
        else
            buffer = buffer - 1;
            writeToMemory = writeToMemory + 1;
            time = time + writeTime(end);
            buffer_with_time = [buffer_with_time; time buffer M_cycles overflow];
        end
        if mod(time, Clk_P) == 0
            break;
        end
    end
    time = ceil(time/Clk_P)*Clk_P;
end

maxBufferSize = max(buffer_with_time(:, 2));

output = {maxBufferSize, buffer_with_time, writeTime};

% fprintf('Mean of Write Times = %f\n', mean(writeTime));
% fprintf('Standard Deviation of Write Times = %f\n', std(writeTime));
% 
% figure
% plot(buffer_with_time(:, 1), buffer_with_time(:, 2))
% xlabel('Time (Clk_{M}s)')
% ylabel('Buffer Size')
% grid on;
% 
% figure
% histogram(writeTime)
% xlabel('Write Time (Clk_{M}s)')
% ylabel('Count')

end