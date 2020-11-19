clc;
close all;
clear;

P_cycles = 20000;

mu = 5;
% muArray = mu_min:2:mu_max;
sigma = 1;
offsetMult_min = 0.5;
offsetMult_max = 2;
offsetMultArray = offsetMult_min:0.2:offsetMult_max;
numIter = 10;
maxBufferArray = zeros(length(offsetMultArray), numIter);
BufferWithTimeArray = {};
writeTimes = {};
leftDiscard = 0;
rightDiscard = 4;


disp('Simulation Running');
for i = 1:1:numIter
    parfor j = 1:1:length(offsetMultArray)
        offsetMult = offsetMultArray(j);
        output = MMU(P_cycles, mu, sigma, leftDiscard, rightDiscard, offsetMult);
        maxBufferArray(j, i) = output{1};
        BufferWithTimeArray{j} = output{2}(:, [1, 2]);
        writeTimes{j} = output{3};
        % disp(j);
    end
    fprintf('%0.1f%% Completed\n', i/numIter*100);
end

fig1 = figure;
errorbar(2./offsetMultArray, mean(maxBufferArray, 2), mean(maxBufferArray, 2) - min(maxBufferArray, [], 2), max(maxBufferArray, [], 2) - mean(maxBufferArray, 2), '-o')
xlabel('offsetMult')
ylabel('Required Buffer Size')
% saveas(fig1, sprintf('Req_BufferSize_with_offsetMult_mu%d_sigma%d.fig', sigma, offsetMult));
% saveas(fig1, sprintf('Req_BufferSize_with_offsetMult_mu%d_sigma%d.png', sigma, offsetMult));

% fig2 = figure;
% for i = 1:1:length(offsetMultArray)
%     plot(BufferWithTimeArray{i}(:, 1), BufferWithTimeArray{i}(:, 2), 'DisplayName', sprintf('offsetMult = %0.1f', offsetMultArray(i)))
%     hold on
% end
% xlabel('Time (Clk_{M}s)')
% ylabel('Buffer Size')
% legend('NumColumns', 2)
% legend show
% grid on;
% saveas(fig2, sprintf('BufferSize_with_time_(different_mu)_sigma%d_offsetMult%0.1f.fig', sigma, offsetMult));
% saveas(fig2, sprintf('BufferSize_with_time_(different_mu)_sigma%d_offsetMult%0.1f.png', sigma, offsetMult));
% 
% fig3 = figure;
% for i = 1:1:length(muArray)
%     histogram(writeTimes{i}, 'DisplayName', sprintf('mu = %d', muArray(i)), 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3)
%     hold on
% end
% xlabel('Time (Clk_{M}s)')
% ylabel('Count')
% legend('NumColumns', 2)
% legend show
% saveas(fig3, sprintf('Histograms_(different_mu)_sigma%d_offsetMult%0.1f.fig', sigma, offsetMult));
% saveas(fig3, sprintf('Histograms_(different_mu)_sigma%d_offsetMult%0.1f.png', sigma, offsetMult));