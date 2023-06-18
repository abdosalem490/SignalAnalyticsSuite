%----------------------------------------%
% Requirement 3
%----------------------------------------%

% creating signal
x = -6:0.01:6;

% quantizing and dequantizing the signal using m = 0
quantizedVals = UniformQuantizer(x, 3, 6, 0);
disp(quantizedVals);
dequantizedVals = UniformDequantizer(quantizedVals, 3, 6, 0);
disp(dequantizedVals);

% plot the graph when m = 0
figure 
title('midrise');
hold on
plot(x, '--','DisplayName','original X value');
plot(dequantizedVals, 'DisplayName','after quantization / dequantization');
legend;
hold off

% quantizing and dequantizing the signal using m = 1
quantizedVals = UniformQuantizer(x, 3, 6, 1);
disp(quantizedVals);
dequantizedVals = UniformDequantizer(quantizedVals, 3, 6, 1);
disp(dequantizedVals);

% plot the graph when m = 1
figure 
title('midtread');
hold on
plot(x, '--','DisplayName','original X value');
plot(dequantizedVals, 'DisplayName','after quantization / dequantization');
legend;
hold off


%----------------------------------------%
% Requirement 4
%----------------------------------------%

% generating random signal
x = -5 + 10 * rand(1, 10000);

% number of bits to be tested on
n_bits = 2:1:8;

% create 2 vectors for storing both theoritcal and actual SNR
SNR_theortical = zeros(1, length(n_bits));
SNR_actual = zeros(1, length(n_bits));

% quantizing/dequantizing over different number of bits but same signal,
% also calculating the theortical 
for i = 1:length(n_bits)
    % quantizing and deqantizing
    quantizedVals = UniformQuantizer(x, n_bits(i), 5, 0);
    dequantizedVals = UniformDequantizer(quantizedVals, n_bits(i), 5, 0);
    % calculate actual SNR
    differenceInSigs = x - dequantizedVals;
    SNR_actual(i) = mean(x.^2) / mean(differenceInSigs.^2); 
    % calculate theoritcal SNR
    LevelsNum = 2^n_bits(i);
    SNR_theortical(i) = 3 * (LevelsNum^2) * mean(x.^2) / (25);
end

% plot the graph , mag2db -> 20log(SNR)
figure 
title('SNR comparison(uniform random input)');
hold on
xlabel('n-bits');
ylabel('SNR (in db)');
plot(n_bits, mag2db(SNR_theortical), '--','DisplayName','theoritical SNR');
plot(n_bits, mag2db(SNR_actual), 'DisplayName','actual SNR');
legend;
hold off

%----------------------------------------%
% Requirement 5
%----------------------------------------%

% array to choose random phase values from 
arr = [-1 1];

% generate the magnitude of the signal
x_mag = exprnd(1, 1, 10000);

% generate phase of the signal
index = randi(2, 1, 10000);
x_phase = arr(index);

% calculate the actual random signal
x = x_mag .* x_phase;

% quantizing/dequantizing over different number of bits but same signal,
% also calculating the SNR theortical 
for i = 1:length(n_bits)
    % quantizing and deqantizing
    quantizedVals = UniformQuantizer(x, n_bits(i), 5, 0);
    dequantizedVals = UniformDequantizer(quantizedVals, n_bits(i), 5, 0);
    % calculate actual SNR
    differenceInSigs = x - dequantizedVals;
    SNR_actual(i) = mean(x.^2) / mean(differenceInSigs.^2); 
    % calculate theoritcal SNR
    LevelsNum = 2^n_bits(i);
    SNR_theortical(i) = 3 * (LevelsNum^2) * mean(x.^2) / (25);
end

% plot the graph , mag2db -> 20log(SNR)
figure 
title('SNR comparison(non-uniform random input)');
hold on
xlabel('n-bits');
ylabel('SNR (in db)');
plot(n_bits, mag2db(SNR_theortical), '--','DisplayName','theoritical SNR');
plot(n_bits, mag2db(SNR_actual), 'DisplayName','actual SNR');
legend;
hold off

%----------------------------------------%
% Requirement 6
%----------------------------------------%

% creating array of possible u
u = [0, 5, 100, 200];

% quantizing/dequantizing over different number of bits but same signal,
% also calculating the SNR theortical 
figure
for i = 1:length(u)
    for j = 1:length(n_bits)

        % get the dequantized signal
        dequantizedVals = compressExpand(x, u(i), n_bits(j), 0);

        % calculate actual SNR
        differenceInSigs = x - dequantizedVals;
        SNR_actual(j) = mean(x.^2) / mean(differenceInSigs.^2);    

        % calculate theoritcal SNR
        LevelsNum = 2^n_bits(j);
        
        if u(i) == 0
            SNR_theortical(j) = 3 * (LevelsNum^2) * mean(x.^2) / (25);
        else
            SNR_theortical(j) = 3 * (LevelsNum^2) / (log(1+u(i)) ^ 2);
        end
        
    end

    % plot the graph , mag2db -> 20log(SNR)
    subplot(2, 2, i); 
    title(['u = ' num2str(u(i))]);
    hold on
    xlabel('n-bits');
    ylabel('SNR (in db)');
    plot(n_bits, mag2db(SNR_theortical), '--','DisplayName','theoritical SNR');
    plot(n_bits, mag2db(SNR_actual), 'DisplayName','actual SNR');
    legend;
    hold off

end


% this fucntion is to compress the signal and return the dequantized value
% after expnading the value again
function [deNormalized] = compressExpand(signal, u, n_bits, m)
    
    % normalize the signal
    normalizedSignal = signal / max(abs(signal));

    % calculate second part of the equation
    compressedVal = (log(1 + u .* abs(normalizedSignal)) / log(1 + u));
    NanIndex = isnan(compressedVal);
    compressedVal(NanIndex) = abs(normalizedSignal(NanIndex));

    % compress the signal
    compressedSignal = (sign(signal) .* compressedVal);
    
    % quantize the signal
    quantizedVals = UniformQuantizer(compressedSignal, n_bits, 1, m);

    % dequantize the signal
    dequantizedVals = UniformDequantizer(quantizedVals, n_bits, 1, m);

    % calculate second part of the equation
    expandedVal = ((((1 + u) .^ abs(dequantizedVals)) - 1) / u);
    NanIndex = isnan(expandedVal);
    expandedVal(NanIndex) = abs(dequantizedVals(NanIndex));

    % expand the signal
    expandedSignal = sign(dequantizedVals) .* expandedVal;
    
    % deNormalize the signal
    deNormalized = expandedSignal .* max(abs(signal));

end


%----------------------------------------%
% Requirement 1
%----------------------------------------%
% this is our quantizer that converts signals to levels (0 -> 2^n-1)
function [q_out] = UniformQuantizer(in_val, n_bits, xmax, m)

    % get the number of the levels
    numOfLevels = 2^n_bits;

    % calculate the step size
    stepSize = (2*xmax) / numOfLevels;
    
    % get the min and max values of the signals
    minVal = (-m * stepSize /2) - xmax + stepSize /2;
    maxVal = (-m * stepSize /2) + xmax - stepSize /2;

    % get the possible the values
    values = minVal:stepSize:maxVal;

    % repeat the vectors so that we can subtract them
    repeatedValues = repmat(values, [length(in_val) 1]);
    repeatedIn_vals = repmat(in_val', [1 length(values)]);

    % subtract the each value in the repeatedIn_val from values to get the
    % least non-negative number 
    subtractedVal = abs(repeatedValues' - repeatedIn_vals');
    
    % get the index of the least non-negative number
    [~, closestIndex] = min(subtractedVal, [], 'omitnan');
    
    % get the closests values to the readings
    q_out = closestIndex - 1;

end



%----------------------------------------%
% Requirement 2
%----------------------------------------%
% this is the dequantizer that converts levels back to amplitude
function [deq_val] = UniformDequantizer(q_ind, n_bits, xmax, m)

    % get the number of the levels
    numOfLevels = 2^n_bits;
    
    % calculate the step size
    stepSize = (2*xmax) / numOfLevels;
 
    % get the min and max values of the signals
    minVal = (-m * stepSize /2) - xmax + stepSize /2;
    
    % calculate the values 
    deq_val = q_ind * stepSize + minVal;
    
end