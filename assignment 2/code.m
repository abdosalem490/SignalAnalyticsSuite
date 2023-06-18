% ------------------------------
% needed variables
% ------------------------------

% this variable is used to switch between plotting the output of each stage
% given N0 = 2 and plotting BER vs E/N0 for different values of N0
displayBlocksOutput = false;

% this variable is used to define the number of bits to be generated
numOfBits = 1000000;

% energy of one symbol = 1
E = 1;

% generate different values for N0 where E/N0  ranges from -10db to 20db
E_over_N0 = -10:0.5:20;   % values in db
N0_arr = E_over_N0 ./ 10;

% select the values of N0 whether it's for plotting blocks output signals
% or plotting BER vs E/N0 for different values of N0
if displayBlocksOutput == true
    N0_arr = 2;
    numOfBits = 10;
else
    N0_arr = E ./ (10 .^ N0_arr);    % getting the actual value of N0
end


simulatedProbabilityOfError1 = zeros(1, length(E_over_N0));  % simulated pobability of error case 1
simulatedProbabilityOfError2 = zeros(1, length(E_over_N0));  % simulated pobability of error case 2
simulatedProbabilityOfError3 = zeros(1, length(E_over_N0));  % simulated pobability of error case 3

theoriticalProbabilityOfError1 = zeros(1, length(E_over_N0)); % theoretical pobability of error case 1
theoriticalProbabilityOfError2 = zeros(1, length(E_over_N0)); % theoretical pobability of error case 2
theoriticalProbabilityOfError3 = zeros(1, length(E_over_N0)); % theoretical pobability of error case 3

% ------------------------------
% BLOCK1 : Binary Source
% ------------------------------

% this is the array that holds our generated pulse
BinarySource = randi([0, 1], 1, numOfBits);

% ------------------------------
% BLOCK2 : pulse shape
% ------------------------------
pulseShapedSignal = pulseShape(BinarySource, 10);

% ------------------------------
% BLOCK3 : channel
% ------------------------------

for i = 1:length(N0_arr)

    % get the current N0
    N0 = N0_arr(i);

    % get the result signal after adding the noise
    R_t = channel(pulseShapedSignal, N0);
    
    % ------------------------------
    % BLOCK4 : Receive filter
    % ------------------------------
    
    % h(t) in case of unit energy filter (case 1)
    h1_t = ones(1, 10);
    convolutedSgn1 = receiveFilter(h1_t, R_t);
    
    % h(t) in case of matched filter doesn't exit (case 2)
    h2_t = [];
    convolutedSgn2 = receiveFilter(h2_t, R_t);
    
    
    % h(t) in case of triangular respose (case 3)
    T = 0:0.11:1;
    h3_t = sqrt(3) .* T;
    convolutedSgn3 = receiveFilter(h3_t, R_t);
    
    % ------------------------------
    % BLOCK5 : sample at T
    % ------------------------------
    
    % this is in case 1 for our matched filter
    sampledData1 = sampler(convolutedSgn1, 10);
    
    
    % this is in case 2 for our matched filter
    sampledData2 = sampler(convolutedSgn2, 10);
    
    
    % this is in case 3 for our matched filter
    sampledData3 = sampler(convolutedSgn3, 10);
    
    
    % ------------------------------
    % BLOCK6 : decode to 0 and 1
    % ------------------------------
    
    
    % this is in case 1 for our matched filter
    decodedVals1 = decode(sampledData1);
    
    % this is in case 2 for our matched filter
    decodedVals2 = decode(sampledData2);
    
    % this is in case 3 for our matched filter
    decodedVals3 = decode(sampledData3);


    % ------------------------------
    % collecting data
    % ------------------------------

    % calculate the simulated probability of error
    simulatedProbabilityOfError1(i) = sum(decodedVals1 ~= BinarySource) / numOfBits;
    simulatedProbabilityOfError2(i) = sum(decodedVals2 ~= BinarySource) / numOfBits;
    simulatedProbabilityOfError3(i) = sum(decodedVals3 ~= BinarySource) / numOfBits;
    
    % importat note: in the attached paper I proved that theoriticalProbabilityOfError2
    % for receiver 2 is equal to 0.5*erfc(1/sqrt(N0)) but the problem in my
    % prove is that I made the area of the impulse signal to be 1, so I
    % made a width of the impulse to be 1/10 of the rectangle so the area
    % of the impulse is 1/10, so after plugining it in the equations, the
    % theoriticalProbabilityOfError2 = 0.5*erfc(1/sqrt(N0*10))

    % calculate the theoretical probability of error
    theoriticalProbabilityOfError1(i) = 0.5*erfc(1/sqrt(N0));
    theoriticalProbabilityOfError2(i) = 0.5*erfc(1/sqrt(N0*10));
    theoriticalProbabilityOfError3(i) = 0.5*erfc(sqrt(3)/(2*sqrt(N0)));

    if displayBlocksOutput ==  true
        % plot the ouput from each stage
        
        % plotting the output of binary source
        T = 0:1:numOfBits-1;    % constructing the time signal
        figure
        title('output of binary source');
        hold on
        xlabel('time');
        ylabel('amplitude');
        h = stem(T, BinarySource, 'DisplayName', 'Binary Source');
        h.Color = [0.5 0 0.8];
        legend
        hold off

        % plotting the output of pulse shape
        T = 0:0.1:(numOfBits-0.1);    % constructing the time signal
        figure
        title('output of pulse shape');
        hold on
        xlabel('time');
        ylabel('amplitude');
        h = stem(T, pulseShapedSignal, 'DisplayName', 'pulse shape');
        h.Color = [0.5 0 0.8];
        legend
        hold off       

        % plotting the output of channel
        T = 0:0.1:(numOfBits-0.1);    % constructing the time signal
        figure
        title('output of the channel');
        hold on
        xlabel('time');
        ylabel('amplitude');
        h = plot(T, R_t, 'DisplayName', 'channel output');
        h.Color = [0.5 0 0.8];
        legend
        hold off       


        % plotting the output of receive filter along with sampling 
        T = 0:0.1:(numOfBits+0.8);    % constructing the time signal
        indexes = NaN(1, numOfBits*11 - 1);
        for i = 10:10:(numOfBits*11 - 1)
            indexes(i) = 1;
        end

        % case 1
        figure
        title('output of the receiver1');
        hold on
        xlabel('time');
        ylabel('amplitude');
        h = plot(T, convolutedSgn1, 'DisplayName', 'h1 output');
        h.Color = [0.5 0 0.8];
        % plotting sampling time
        sampleFunc = indexes .* convolutedSgn1;
        stem(T, sampleFunc, 'DisplayName', 'sample at T' );
        legend
        hold off 
        
        % case 2
        figure
        title('output of the receiver2');
        hold on
        xlabel('time');
        ylabel('amplitude');
        h = plot(T, convolutedSgn2, 'DisplayName', 'h2 output');
        h.Color = [0.5 0 0.8];
        % plotting sampling time
        sampleFunc = indexes .* convolutedSgn2;
        stem(T, sampleFunc, 'DisplayName', 'sample at T' );
        legend
        hold off 

        % case 3
        figure
        title('output of the receiver3');
        hold on
        xlabel('time');
        ylabel('amplitude');
        h = plot(T, convolutedSgn3, 'DisplayName', 'h3 output');
        h.Color = [0.5 0 0.8];
        % plotting sampling time
        sampleFunc = indexes .* convolutedSgn3;
        stem(T, sampleFunc, 'DisplayName', 'sample at T' );
        legend
        hold off

        % plotting the output of decoders
        T = 0:1:numOfBits-1;    % constructing the time signal

        % case 1
        figure
        title('output of the decoder1');
        hold on
        xlabel('time');
        ylabel('amplitude');
        h = stem(T, decodedVals1, 'DisplayName', 'decoder1 output');
        h.Color = [0.5 0 0.8];
        legend
        hold off       

        % case 2
        figure
        title('output of the decoder2');
        hold on
        xlabel('time');
        ylabel('amplitude');
        h = stem(T, decodedVals2, 'DisplayName', 'decoder2 output');
        h.Color = [0.5 0 0.8];
        legend
        hold off

        % case 3
        figure
        title('output of the decoder3');
        hold on
        xlabel('time');
        ylabel('amplitude');
        h = stem(T, decodedVals3, 'DisplayName', 'decoder3 output');
        h.Color = [0.5 0 0.8];
        legend
        hold off
        
    end
    
end 

% ------------------------------
% plotting BER vs E/N0 for different values of N0
% ------------------------------
if displayBlocksOutput == false
    
    % plotting the theoretical and simulated Bit Error Rate (BER) Vs
    % E/No for receiver 1
    figure
    title('BER vs E/N for receiver1');
    hold on
    xlabel('E/N0');
    ylabel('BER');
    set(gca, 'YScale', 'log')
    h = plot(E_over_N0, simulatedProbabilityOfError1, 'DisplayName', 'simulated BER');
    h.Color = [0.5 0 0.8];
    plot(E_over_N0, theoriticalProbabilityOfError1, '--', 'DisplayName', 'theoretical BER');
    legend
    hold off

    % plotting the theoretical and simulated Bit Error Rate (BER) Vs
    % E/No for receiver 2
    figure
    title('BER vs E/N for receiver2');
    hold on
    xlabel('E/N0');
    ylabel('BER');
    set(gca, 'YScale', 'log')
    h = plot(E_over_N0, simulatedProbabilityOfError2, 'DisplayName', 'simulated BER');
    h.Color = [0.5 0 0.8];
    plot(E_over_N0, theoriticalProbabilityOfError2, '--', 'DisplayName', 'theoretical BER');
    legend
    hold off

    % plotting the theoretical and simulated Bit Error Rate (BER) Vs
    % E/No for receiver 3
    figure
    title('BER vs E/N for receiver3');
    hold on
    xlabel('E/N0');
    ylabel('BER');
    set(gca, 'YScale', 'log')
    h = plot(E_over_N0, simulatedProbabilityOfError3, 'DisplayName', 'simulated BER');
    h.Color = [0.5 0 0.8];
    plot(E_over_N0, theoriticalProbabilityOfError3, '--', 'DisplayName', 'theoretical BER');
    legend
    hold off


end    



% this function is to convert pulses to (+ and -) rectangle pulses (pulse
% shape blocl)
function [G_t] = pulseShape(binaryBitsArr, numOfPulsesPerOnePulse)
    
    % needed varaibles for our for loop
    T = 0:1:numOfPulsesPerOnePulse-1;
    T = T / numOfPulsesPerOnePulse;
    V = [];

    % repeat each pulse 10 times to make a rectangle shape pulse
    for bit = binaryBitsArr
        % convert each signal to +1 and -1 insteal of 0 and 1
        if bit == 0
            V = cat(2, V, -rectpuls(T, numOfPulsesPerOnePulse));
        else
            V = cat(2, V, rectpuls(T, numOfPulsesPerOnePulse));
        end
    end

    % return the result
    G_t = V;
end


% this function is to add noise to the channel (N0)
function [R_t] = channel(originalSignal, N0)
    
    % generate the noise from gaussian distribution given mean = 0 and
    % varience = N0 / 2, refer to https://ch.mathworks.com/help/matlab/math/random-numbers-with-specific-mean-and-variance.html  
    noise = sqrt(N0 / 2) .* randn([1, length(originalSignal)]);

    % add the noise to the original signal (sqrt(10) is to normalized the length of the one signal pulse)
    R_t = originalSignal/sqrt(10) + noise;

end

% this function is used to convolute the matched filter with the signal
% output from the channel
function [convolutedSgn] = receiveFilter(H_t, R_t) 

    if(isempty(H_t))
        % the matched filter doesn't exist
        convolutedSgn = cat(2, R_t, zeros(1, 9));
    else
        % convulte the signal output from the channel with matched filter
        % response
        convolutedSgn = conv(H_t, R_t);
    end

end


% this function is used to sample at time T (sample every N values)
function [sampledData] = sampler(data, N)    
    
    % temporary vector
    tempVec = zeros(1, floor(length(data)/N));

    % sample every N times
    for i = 10:N:length(data)
        tempVec(i/10) = data(i);
    end
    
    % return the result
    sampledData = tempVec;

end


% this function is used to decode the results to either 0 or 1 
function [decodedVals] = decode(data)
    % create temporary vector
    tempVec = data;

    % map the values to either 1 or 0
    tempVec(tempVec <= 0) = 0;
    tempVec(tempVec > 0) = 1;

    % return the result
    decodedVals = tempVec;
end