%% 
clc; clear;

modulationSchemes = [4, 16, 64, 256];
numInputBits = input('Enter total number of input bits: ');
SNR_range_dB = 0:0.5:30;

berMatrix_equal = zeros(length(modulationSchemes), length(SNR_range_dB));
berMatrix_water = zeros(length(modulationSchemes), length(SNR_range_dB));

for m = 1:length(modulationSchemes)
    for useWaterFilling = [0, 1]
        modOrder = modulationSchemes(m);
        fprintf('\n===== Running for %d-QAM with%s Water-Filling =====\n', modOrder, ternary(useWaterFilling, '', 'out'));

        %% Parameters
        subcarrierSpacing = 30e3;
        bandwidth = 100e6;
        numSlots = 2;
        N_fft = 4096;
        numUsedSubcarriers = 3300;
        cpDuration_us = 2.3;
        addNoise = true;

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

        berVals = zeros(size(SNR_range_dB));

        for snrIdx = 1:length(SNR_range_dB)
            SNR_dB = SNR_range_dB(snrIdx);
            dataBits = randi([0, 1], totalBitsUsed, 1);
            qamSymbols = qammod(dataBits, modOrder, 'InputType', 'bit', 'UnitAveragePower', true);
            parallelData = reshape(qamSymbols, numUsedSubcarriers, totalSymbols);

            subcarrierSNR_dB = 5 + 25 * rand(numUsedSubcarriers, 1);
            subcarrierNoise = 10.^(-subcarrierSNR_dB / 10);
            totalPower = numUsedSubcarriers;

            objectiveFn = @(level) abs(sum(max(level - subcarrierNoise, 0)) - totalPower);
            lowerBound = min(subcarrierNoise);
            upperBound = max(subcarrierNoise) + totalPower;
            waterLevel = fminbnd(objectiveFn, lowerBound, upperBound);
            allocatedPower = max(waterLevel - subcarrierNoise, 0);

            if useWaterFilling
                powerScaling = sqrt(allocatedPower);
            else
                powerScaling = ones(size(allocatedPower));
            end

            scaledSymbols = parallelData .* powerScaling;

            ofdmGrid = zeros(N_fft, totalSymbols);
            ofdmGrid(centerStart:centerStart+numUsedSubcarriers-1, :) = scaledSymbols;
            ifftData = manual_radix2_ifft(ifftshift(ofdmGrid,1));
            serialData = ifftData(:);

            reshaped = reshape(serialData, N_fft, totalSymbols);
            ofdmSignal = [];
            for i = 1:totalSymbols
                frame = reshaped(:,i);
                cp = frame(end-cpLength+1:end);
                ofdmSignal = [ofdmSignal; cp; frame];
            end

            if addNoise
                signalPower = mean(abs(ofdmSignal).^2);
                noisePower = signalPower / (10^(SNR_dB/10));
                noise = sqrt(noisePower/2) * (randn(size(ofdmSignal)) + 1j*randn(size(ofdmSignal)));
                rx_serial = ofdmSignal + noise;
            else
                rx_serial = ofdmSignal;
            end

            rx_matrix = reshape(rx_serial, N_fft + cpLength, totalSymbols);
            rx_noCP = rx_matrix(cpLength+1:end, :);
            rx_fft = fftshift(fft(rx_noCP, N_fft, 1), 1);
            rx_subcarriers = rx_fft(centerStart:centerStart+numUsedSubcarriers-1, :);
            rx_serial_symbols = rx_subcarriers(:);
            rx_eq = rx_serial_symbols ./ repmat(powerScaling, totalSymbols, 1);
            rx_bits = qamdemod(rx_eq, modOrder, 'OutputType', 'bit', 'UnitAveragePower', true);

            effectiveBits = min(numInputBits, length(rx_bits));
            bitErrors = sum(rx_bits(1:effectiveBits) ~= dataBits(1:effectiveBits));
            berVals(snrIdx) = bitErrors / numInputBits;

            %fprintf('Modulation: %d-QAM | SNR = %.1f dB -> BER = %.5e\n', modOrder, SNR_dB, berVals(snrIdx));

            if m == length(modulationSchemes) && snrIdx == length(SNR_range_dB)
                dataBits_store = dataBits;
                rxBits_store = rx_bits;
                ofdmSignal_store = ofdmSignal;
                rxSerial_store = rx_serial;
                allocatedPower_store = allocatedPower;
            end
        end

        if useWaterFilling
            berMatrix_water(m, :) = berVals;
        else
            berMatrix_equal(m, :) = berVals;
        end
    end
end

% Time vector for plotting
timeVec = (0:length(ofdmSignal_store)-1) * Ts * 1e6;

% Original Plots for last modulation
figure;
subplot(4,1,1);
plot(timeVec, real(ofdmSignal_store)); hold on;
plot(timeVec, imag(ofdmSignal_store));
title('Transmitted OFDM Signal (Time Domain)');
xlabel('Time [\mus]'); ylabel('Amplitude'); legend('Real','Imag');

subplot(4,1,2);
plot(timeVec, real(rxSerial_store)); hold on;
plot(timeVec, imag(rxSerial_store));
title('Received OFDM Signal (Time Domain)');
xlabel('Time [\mus]'); ylabel('Amplitude'); legend('Real','Imag');

subplot(4,1,3);
stem(dataBits_store(1:100), 'filled');
title('Original Input Bits (First 100)');
xlabel('Bit Index'); ylabel('Bit Value');

subplot(4,1,4);
stem(rxBits_store(1:100), 'filled');
title('Received Output Bits (First 100)');
xlabel('Bit Index'); ylabel('Bit Value');

figure;
subplot(2,1,1);
stem(allocatedPower_store, 'filled');
title('Power Allocated to Subcarriers (Water-Filling)');
xlabel('Subcarrier Index'); ylabel('Power');

subplot(2,1,2);
semilogy(SNR_range_dB, berMatrix_water(end, :), 'o-', 'LineWidth', 1.5);
grid on;
title(sprintf('BER vs. SNR for %d-QAM (5G NR-style OFDM)', modulationSchemes(end)));
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');

figure;
colors = ['r', 'g', 'b', 'k'];
for m = 1:length(modulationSchemes)
    semilogy(SNR_range_dB, berMatrix_equal(m, :), ['--' colors(m)], 'LineWidth', 1.2); hold on;
    semilogy(SNR_range_dB, berMatrix_water(m, :), ['-' colors(m)], 'LineWidth', 1.5);
end
legend({'QPSK (Equal)', 'QPSK (Water)', '16-QAM (Equal)', '16-QAM (Water)', '64-QAM (Equal)', '64-QAM (Water)', '256-QAM (Equal)', '256-QAM (Water)'}, 'Location', 'southwest');
grid on;
title('BER vs. SNR for Multiple Modulation Schemes (Equal vs Water-Filling)');
xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');

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

function out = ternary(cond, valTrue, valFalse)
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end
