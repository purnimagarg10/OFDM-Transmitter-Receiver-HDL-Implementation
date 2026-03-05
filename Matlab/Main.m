% ======================== MAIN SCRIPT ========================
clc; clear;

%% ------------------------- User Inputs -------------------------
modInput = input('Enter modulation scheme (QPSK=4, 16, 64, 256): ');
if ~ismember(modInput, [4, 16, 64, 256])
    error('Invalid modulation scheme selected.');
end
modOrder = modInput;

numInputBits = input('Enter total number of input bits: ');
useWaterFilling = input('Use water-filling? (1 = Yes, 0 = No): ');
addAWGN = input('Add AWGN channel? (1 = Yes, 0 = No): ');
useTDL = input('Use TDL channel model? (1 = Yes, 0 = No): ');

%% ------------------------- Parameters -------------------------
subcarrierSpacing = 30e3;
bandwidth = 100e6;
numSlots = 2;
N_fft = 4096;
numUsedSubcarriers = 3300;
cpDuration_us = 2.3;
SNR_range_dB = 0:0.5:40;

numTx = 2;
numRx = 2;

%% ------------------------- Derived -------------------------
bitsPerSymbol = log2(modOrder);
symbolsPerSlot = ceil(numInputBits / (bitsPerSymbol * numUsedSubcarriers));
totalSymbols = symbolsPerSlot * numSlots;
totalBitsUsed = totalSymbols * numUsedSubcarriers * bitsPerSymbol;

T_sym = 1 / subcarrierSpacing;
fs = subcarrierSpacing * N_fft;
Ts = 1 / fs;
cpRatio = (cpDuration_us * 1e-6) / T_sym;
cpLength = round(cpRatio * N_fft);
centerStart = floor(N_fft/2) - floor(numUsedSubcarriers/2) + 1;

fprintf('\nSelected Modulation: %d-QAM\n', modOrder);
fprintf('Input Bits: %d\n', numInputBits);
fprintf('Sampling Frequency = %.2f MHz\n', fs/1e6);
fprintf('Sampling Time = %.3f ns\n', Ts*1e9);
fprintf('Cyclic Prefix Length = %d samples\n\n', cpLength);

berVals = zeros(size(SNR_range_dB));

numTrials = 10;  % Number of channel realizations per SNR point

for snrIdx = 1:length(SNR_range_dB)
    SNR_dB = SNR_range_dB(snrIdx);
    snrLinear = 10^(SNR_dB/10);

    totalBitErrors = 0;
    totalBitsEvaluated = 0;

    for trial = 1:numTrials
        dataBits = generate_random_bits(totalBitsUsed);
        qamSymbols = qam_modulation(dataBits, modOrder);
        parallelData = reshape(qamSymbols, numUsedSubcarriers, totalSymbols);

        if useWaterFilling
            subcarrierSNR_dB = 5 + 45 * rand(numUsedSubcarriers, 1);
            subcarrierNoise = 10.^(-subcarrierSNR_dB / 10);
        else
            subcarrierNoise = ones(numUsedSubcarriers, 1);
        end

        totalPower = numUsedSubcarriers;
        powerScaling = apply_water_filling(subcarrierNoise, totalPower, useWaterFilling);
        scaledSymbols = parallelData .* powerScaling;

        ofdmGrid = zeros(N_fft, totalSymbols);
        ofdmGrid(centerStart:centerStart+numUsedSubcarriers-1, :) = scaledSymbols;
        ifftData = manual_radix2_ifft(ifftshift(ofdmGrid, 1));

        ofdmWithCP = add_guard_interval(ifftData, cpLength, N_fft, totalSymbols);
        txSignalMIMO = repmat(ofdmWithCP, 1, numTx).';

        if useTDL
            [rxMIMO, hFreq] = tdl_channel_model_fd(txSignalMIMO, N_fft, cpLength, numUsedSubcarriers, totalSymbols, numRx, numTx);
            mimoChannelFD = hFreq;
        else
            if addAWGN
                mimoChannel = (randn(numRx, numTx) + 1j * randn(numRx, numTx)) / sqrt(2);
            else
                mimoChannel = eye(numRx, numTx);
            end
            rxMIMO = mimoChannel * txSignalMIMO;
            rxMIMO = apply_channel(rxMIMO, SNR_dB, addAWGN);
        end

        rx_eq = zeros(numUsedSubcarriers, totalSymbols);
        for i = 1:totalSymbols
            symbolMat = zeros(numRx, N_fft + cpLength);
            for r = 1:numRx
                startIdx = (i-1)*(N_fft + cpLength) + 1;
                endIdx = i*(N_fft + cpLength);
                symbolMat(r,:) = rxMIMO(r, startIdx:endIdx);
            end

            symbolMat = symbolMat(:, cpLength+1:end);
            fftData = fftshift(fft(symbolMat, N_fft, 2), 2);
            rxSub = fftData(:, centerStart:centerStart+numUsedSubcarriers-1);

            for sc = 1:numUsedSubcarriers
                y = rxSub(:, sc);
                if useTDL
                    Hf = squeeze(mimoChannelFD(:,:,sc));
                    mmseEq = (Hf' * Hf + (1/snrLinear) * eye(numTx)) \ Hf';
                elseif addAWGN
                    mmseEq = (mimoChannel' * mimoChannel + (1/snrLinear) * eye(numTx)) \ mimoChannel';
                else
                    mmseEq = pinv(mimoChannel);
                end
                rx_eq(sc, i) = mmseEq(1,:) * y;
            end
        end

        rx_eq = rx_eq(:) ./ repmat(powerScaling, totalSymbols, 1);
        rx_bits = qamdemod(rx_eq, modOrder, 'OutputType', 'bit', 'UnitAveragePower', true);

        effectiveBits = min(numInputBits, length(rx_bits));
        bitErrors = sum(rx_bits(1:effectiveBits) ~= dataBits(1:effectiveBits));

        totalBitErrors = totalBitErrors + bitErrors;
        totalBitsEvaluated = totalBitsEvaluated + effectiveBits;

        if trial == 1 && snrIdx == length(SNR_range_dB)
            dataBits_store = dataBits;
            rxBits_store = rx_bits;
            power_store = powerScaling.^2;
            ofdmSignal_store = ofdmWithCP;
            rxSerial_store = rxMIMO(1,:).';
            timeVec = (0:length(ofdmWithCP)-1) * Ts * 1e6;
        end
    end

    berVals(snrIdx) = totalBitErrors / totalBitsEvaluated;
    fprintf('SNR = %2.1f dB -> Avg BER = %.6f\n', SNR_dB, berVals(snrIdx));
end


%% ------------------------- Plots -------------------------
figure;
subplot(3,1,1); stem(dataBits_store(1:100), 'filled');
title('Original Input Bits'); xlabel('Bit Index'); ylabel('Bit Value');

subplot(3,1,2); stem(rxBits_store(1:100), 'filled');
title('Received Output Bits'); xlabel('Bit Index'); ylabel('Bit Value');

subplot(3,1,3); semilogy(SNR_range_dB, berVals, 'o-', 'LineWidth', 1.5); grid on;
title(sprintf('BER vs. SNR (%d-QAM, MIMO)', modOrder)); xlabel('SNR (dB)'); ylabel('BER');

figure;
subplot(3,1,1);
plot(real(ofdmSignal_store)); hold on; plot(imag(ofdmSignal_store));
title('Transmitted OFDM Signal (Samples)'); legend('Real','Imag'); ylabel('Amplitude');

subplot(3,1,2);
plot(real(rxSerial_store)); hold on; plot(imag(rxSerial_store));
title('Received OFDM Signal (Samples)'); legend('Real','Imag'); ylabel('Amplitude');

subplot(3,1,3);
stem(power_store, 'filled'); title('Power Allocation per Subcarrier');
xlabel('Subcarrier Index'); ylabel('Power');

figure;
subplot(2,1,1);
plot(timeVec, real(ofdmSignal_store)); hold on;
plot(timeVec, imag(ofdmSignal_store));
title('Transmitted OFDM Signal vs Time'); xlabel('Time [\mus]'); ylabel('Amplitude'); legend('Real','Imag');

subplot(2,1,2);
plot(timeVec, real(rxSerial_store)); hold on;
plot(timeVec, imag(rxSerial_store));
title('Received OFDM Signal vs Time'); xlabel('Time [\mus]'); ylabel('Amplitude'); legend('Real','Imag');

%% ------------------------- Supporting Functions -------------------------
function bits = generate_random_bits(numBits)
    bits = randi([0, 1], numBits, 1);
end

function symbols = qam_modulation(bits, modOrder)
    symbols = qammod(bits, modOrder, 'InputType', 'bit', 'UnitAveragePower', true);
end

function scaling = apply_water_filling(noise, totalPower, useWF)
    if useWF
        objectiveFn = @(level) abs(sum(max(level - noise, 0)) - totalPower);
        lowerBound = min(noise);
        upperBound = max(noise) + totalPower;
        waterLevel = fminbnd(objectiveFn, lowerBound, upperBound);
        allocatedPower = max(waterLevel - noise, 0);
        scaling = sqrt(allocatedPower);
    else
        scaling = ones(size(noise));
    end
end

function rx = apply_channel(tx, SNR_dB, addNoise)
    if addNoise
        signalPower = mean(abs(tx(:)).^2);
        noisePower = signalPower / (10^(SNR_dB/10));
        noise = sqrt(noisePower/2) * (randn(size(tx)) + 1j*randn(size(tx)));
        rx = tx + noise;
    else
        rx = tx;
    end
end

function x = manual_radix2_ifft(X)
    N = size(X, 1);
    X = bit_reverse_order(X, N);
    stages = log2(N);
    x = X;
    for stage = 1:stages
        numButterflies = 2^(stage-1);
        butterflySeparation = 2^stage;
        W = exp(2j * pi * (0:numButterflies-1) / butterflySeparation);
        for k = 0:(N/butterflySeparation - 1)
            for j = 0:(numButterflies - 1)
                idx1 = k * butterflySeparation + j + 1;
                idx2 = idx1 + numButterflies;
                temp = x(idx2,:) .* W(j+1);
                x(idx2,:) = x(idx1,:) - temp;
                x(idx1,:) = x(idx1,:) + temp;
            end
        end
    end
    x = x / N;
end

function Xr = bit_reverse_order(X, N)
    Xr = zeros(size(X));
    numBits = log2(N);
    for i = 0:N-1
        revIdx = bin2dec(fliplr(dec2bin(i, numBits))) + 1;
        Xr(revIdx,:) = X(i+1,:);
    end
end

function ofdmSignal = add_guard_interval(serialData, cpLength, N_fft, totalSymbols)
    reshaped = reshape(serialData, N_fft, totalSymbols);
    ofdmSignal = [];
    for i = 1:totalSymbols
        frame = reshaped(:,i);
        cp = frame(end-cpLength+1:end);
        ofdmSignal = [ofdmSignal; cp; frame];
    end
end

