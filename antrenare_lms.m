clc
clear
close all

%% ANALIZĂ SUPLIMENTARĂ:
% Influența numărului de simboluri asupra performanței metodei LMS

%% PARAMETRI SISTEM

M = 16;                  % Modulație 16-QAM
k = log2(M);             % Biți per simbol
Fs = 10e3;               % Rată de eșantionare [Hz]
SNR_dB = 15;             % Raport semnal-zgomot [dB]

% Parametri canal
fadingStrength = 0.1;    % Intensitate fading
phaseOffset    = 0;      % Offset de fază [rad]

% Oversampling și pulse shaping
sps     = 4;             % Samples per symbol
rolloff = 0.25;          % Raised cosine roll-off
span    = 6;             % Filter span în simboluri

% Parametri LMS
mu = 0.01;
N_lms = 11;

% Pentru comparație corectă, training-ul este păstrat proporțional cu N
trainingRatio = 0.02;    % 2% din simboluri pentru antrenare

% Valorile analizate pentru numărul de simboluri
N_values = [10000 25000 50000 100000 250000];

% Etichete pentru axa OX, fără notația ×10^5
N_labels = {'10.000','25.000','50.000','100.000','250.000'};

% Filtru SRRC
rcFilt = rcosdesign(rolloff, span, sps, 'sqrt');

% Alocare rezultate
EVM_LMS_N = zeros(length(N_values),1);
BER_LMS_N = zeros(length(N_values),1);
SNR_LMS_N = zeros(length(N_values),1);

N_training_values = zeros(length(N_values),1);
N_bits_values = zeros(length(N_values),1);

maxShift = 10;

%% RULARE AUTOMATĂ PENTRU FIECARE VALOARE N

fprintf('\n============================================================\n');
fprintf(' ANALIZA LMS ÎN FUNCȚIE DE NUMĂRUL DE SIMBOLURI\n');
fprintf('============================================================\n');
fprintf('N simboluri\tN biți\t\tN training\tEVM(%%)\t\tBER\t\tSNR(dB)\n');

for idxN = 1:length(N_values)

    N = N_values(idxN);
    N_bits_values(idxN) = N * k;

    % Reproducibilitate pentru fiecare rulare
    rng(42);

    %% GENERARE SEMNAL TX

    data = randi([0 1], N*k, 1);

    txSig = qammod(data, M, ...
        'InputType','bit', ...
        'UnitAveragePower',true);

    %% OVERSAMPLING + PULSE SHAPING

    txUp = upsample(txSig, sps);
    txWave = conv(txUp, rcFilt, 'same');

    %% CANAL: FADING + AWGN

    fadingCoeff_wave = (randn(size(txWave)) + 1j*randn(size(txWave))) / sqrt(2);
    fadingCoeff_wave = (1 - fadingStrength) + fadingStrength * fadingCoeff_wave;

    rxWave = txWave .* fadingCoeff_wave * exp(1j*phaseOffset);
    rxWave = awgn(rxWave, SNR_dB, 'measured');

    %% LMS: MATCHED FILTER + DOWNSAMPLING

    rxMatched_LMS = conv(rxWave, rcFilt, 'same');

    delay_LMS = span*sps/2;
    rxDownsampled_LMS = rxMatched_LMS(delay_LMS+1:sps:end);

    Nsym_LMS = min(length(rxDownsampled_LMS), N);

    rxDownsampled_LMS = rxDownsampled_LMS(1:Nsym_LMS);
    txAligned_LMS = txSig(1:Nsym_LMS);

    %% COMPENSARE INIȚIALĂ DE FAZĂ

    phaseEst_LMS_init = angle(mean((rxDownsampled_LMS.^4) .* ...
                            conj(txAligned_LMS.^4))) / 4;

    rxDownsampled_LMS = rxDownsampled_LMS * exp(-1j*phaseEst_LMS_init);

    %% NORMALIZARE

    rxDownsampled_LMS = rxDownsampled_LMS / rms(rxDownsampled_LMS);

    %% PARAMETRI TRAINING PENTRU N CURENT

    N_training = round(trainingRatio * Nsym_LMS);
    N_training = max(N_training, N_lms + 1);

    N_training_values(idxN) = N_training;

    %% ALGORITM LMS

    w = zeros(N_lms, 1);
    w(ceil(N_lms/2)) = 1;

    y_lms = zeros(Nsym_LMS, 1);

    rxPadded = [zeros(N_lms-1, 1); rxDownsampled_LMS];

    % Faza de antrenare
    for n = 1:N_training

        idx = n + N_lms - 1;
        x_vec = rxPadded(idx:-1:idx-N_lms+1);

        y_lms(n) = w' * x_vec;

        e = txAligned_LMS(n) - y_lms(n);

        w = w + mu * conj(e) * x_vec;
    end

    % Regim decision-directed
    for n = N_training+1:Nsym_LMS

        idx = n + N_lms - 1;
        x_vec = rxPadded(idx:-1:idx-N_lms+1);

        y_lms(n) = w' * x_vec;

        decision = qammod(qamdemod(y_lms(n), M, ...
                          'UnitAveragePower', true), ...
                          M, 'UnitAveragePower', true);

        e = decision - y_lms(n);

        w = w + mu * conj(e) * x_vec;
    end

    %% NORMALIZARE FINALĂ

    y_lms = y_lms / rms(y_lms);

    %% CALCUL METRICI

    [EVM_LMS_N(idxN), BER_LMS_N(idxN), SNR_LMS_N(idxN)] = ...
        metricsQAM_local(txSig, y_lms, data, M, k, maxShift);

    %% AFIȘARE ÎN COMMAND WINDOW

    fprintf('%d\t\t%d\t\t%d\t\t%.2f\t\t%.6f\t%.2f\n', ...
        N, N_bits_values(idxN), N_training, ...
        EVM_LMS_N(idxN), BER_LMS_N(idxN), SNR_LMS_N(idxN));

end

fprintf('============================================================\n');

%% TABEL FINAL

Rezultate_LMS_N = table( ...
    N_values(:), ...
    N_bits_values(:), ...
    N_training_values(:), ...
    EVM_LMS_N, ...
    BER_LMS_N, ...
    SNR_LMS_N, ...
    'VariableNames', {'N_simboluri', 'N_biti', 'N_training', ...
                      'EVM_procente', 'BER', 'SNR_dB'});

disp(' ');
disp('TABEL FINAL - Influența numărului de simboluri asupra LMS');
disp(Rezultate_LMS_N);

% Export în Excel
writetable(Rezultate_LMS_N, 'Rezultate_LMS_in_functie_de_N.xlsx');

%% FIGURA 1: EVM LMS ÎN FUNCȚIE DE N

figure('Name','LMS - EVM în funcție de N', ...
       'Position',[300 200 850 500]);

plot(N_values, EVM_LMS_N, '-o', ...
     'LineWidth',1.8, ...
     'MarkerSize',7);

xlabel('Număr de simboluri N');
ylabel('EVM (%)');
title('Influența numărului de simboluri asupra EVM pentru metoda LMS', ...
      'FontWeight','bold');

grid on;
set(gca,'FontSize',10);

xticks(N_values);
xticklabels(N_labels);
ax = gca;
ax.XAxis.Exponent = 0;

xlim([min(N_values)*0.9 max(N_values)*1.05]);

for i = 1:length(N_values)
    text(N_values(i), EVM_LMS_N(i), sprintf('%.2f%%', EVM_LMS_N(i)), ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','center', ...
        'FontSize',9);
end

%% FIGURA 2: BER LMS ÎN FUNCȚIE DE N

figure('Name','LMS - BER în funcție de N', ...
       'Position',[350 220 850 500]);

BER_LMS_plot = BER_LMS_N;

% Dacă BER este zero, pe scară logaritmică se afișează limita statistică 1/(N*k)
for i = 1:length(BER_LMS_plot)
    if BER_LMS_plot(i) == 0
        BER_LMS_plot(i) = 1 / (N_values(i) * k);
    end
end

semilogy(N_values, BER_LMS_plot, '-o', ...
         'LineWidth',1.8, ...
         'MarkerSize',7);

xlabel('Număr de simboluri N');
ylabel('BER');
title('Influența numărului de simboluri asupra BER pentru metoda LMS', ...
      'FontWeight','bold');

grid on;
set(gca,'FontSize',10);

xticks(N_values);
xticklabels(N_labels);
ax = gca;
ax.XAxis.Exponent = 0;

xlim([min(N_values)*0.9 max(N_values)*1.05]);

for i = 1:length(N_values)

    if BER_LMS_N(i) == 0
        labelBER = sprintf('< %.1e', 1/(N_values(i)*k));
    else
        labelBER = sprintf('%.2e', BER_LMS_N(i));
    end

    text(N_values(i), BER_LMS_plot(i)*1.15, labelBER, ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','center', ...
        'FontSize',9);
end

%% FIGURA 3: SNR LMS ÎN FUNCȚIE DE N

figure('Name','LMS - SNR în funcție de N', ...
       'Position',[400 240 850 500]);

plot(N_values, SNR_LMS_N, '-o', ...
     'LineWidth',1.8, ...
     'MarkerSize',7);

xlabel('Număr de simboluri N');
ylabel('SNR (dB)');
title('Influența numărului de simboluri asupra SNR pentru metoda LMS', ...
      'FontWeight','bold');

grid on;
set(gca,'FontSize',10);

xticks(N_values);
xticklabels(N_labels);
ax = gca;
ax.XAxis.Exponent = 0;

xlim([min(N_values)*0.9 max(N_values)*1.05]);

for i = 1:length(N_values)
    text(N_values(i), SNR_LMS_N(i), sprintf('%.2f dB', SNR_LMS_N(i)), ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','center', ...
        'FontSize',9);
end

%% FUNCȚIE LOCALĂ PENTRU CALCULUL METRICILOR

function [bestEVM, bestBER, bestSNR] = metricsQAM_local(txSig, rxSig, dataBitsAll, M, k, maxShift)

    bestEVM = inf;
    bestBER = NaN;
    bestSNR = NaN;

    txSig = txSig(:);
    rxSig = rxSig(:);
    dataBitsAll = dataBitsAll(:);

    for shift = -maxShift:maxShift

        if shift >= 0
            txStart = 1 + shift;
            rxStart = 1;
        else
            txStart = 1;
            rxStart = 1 - shift;
        end

        L = min(length(txSig) - txStart + 1, length(rxSig) - rxStart + 1);

        if L <= 0
            continue;
        end

        tx = txSig(txStart:txStart+L-1);
        rx = rxSig(rxStart:rxStart+L-1);

        % Corecție complexă amplitudine + fază
        alpha = (rx' * tx) / (rx' * rx + eps);
        rxCorr = alpha * rx;

        err = rxCorr - tx;

        EVM = 100 * rms(err) / rms(tx);

        Ps = mean(abs(tx).^2);
        Pn = mean(abs(err).^2);
        SNR = 10 * log10(Ps / (Pn + eps));

        % Demodulare
        rxBits = qamdemod(rxCorr, M, ...
            'OutputType','bit', ...
            'UnitAveragePower',true);

        rxBits = rxBits(:);

        bitStart = (txStart - 1) * k + 1;
        bitEnd = bitStart + length(rxBits) - 1;

        if bitEnd > length(dataBitsAll)
            Lbits = length(dataBitsAll) - bitStart + 1;
            rxBits = rxBits(1:Lbits);
            txBits = dataBitsAll(bitStart:end);
        else
            txBits = dataBitsAll(bitStart:bitEnd);
        end

        BER = mean(rxBits ~= txBits);

        % Se păstrează alinierea cu EVM minim
        if EVM < bestEVM
            bestEVM = EVM;
            bestBER = BER;
            bestSNR = SNR;
        end
    end

    % Caz de rezervă
    if isinf(bestEVM)

        L = min(length(txSig), length(rxSig));
        tx = txSig(1:L);
        rx = rxSig(1:L);

        alpha = (rx' * tx) / (rx' * rx + eps);
        rxCorr = alpha * rx;

        err = rxCorr - tx;

        bestEVM = 100 * rms(err) / rms(tx);

        Ps = mean(abs(tx).^2);
        Pn = mean(abs(err).^2);
        bestSNR = 10 * log10(Ps / (Pn + eps));

        rxBits = qamdemod(rxCorr, M, ...
            'OutputType','bit', ...
            'UnitAveragePower',true);

        rxBits = rxBits(:);

        Lbits = min(length(rxBits), length(dataBitsAll));
        bestBER = mean(rxBits(1:Lbits) ~= dataBitsAll(1:Lbits));
    end
end