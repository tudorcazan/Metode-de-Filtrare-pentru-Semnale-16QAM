clc; clear; close all;

%% PARAMETRI SISTEM
M = 16;
N = 250000;
k = log2(M);
Fs = 10e3;

% Parametri canal
fadingStrength = 0;
phaseOffset = 0;

% Parametri oversampling și pulse shaping
sps = 4;
rolloff = 0.25;
span = 6;

%% PARAMETRI SIMULARE SNR
SNR_range = 0:1:20;
num_SNR = length(SNR_range);

% Praguri pentru evaluarea performanței
BER_threshold = 1e-4;
BER_threshold_excellent = 1e-6;
EVM_threshold = 12.5;
SNR_threshold = 15;

%% INIȚIALIZARE MATRICE REZULTATE
EVM_results = zeros(num_SNR, 5);
BER_results = zeros(num_SNR, 5);
SNR_results = zeros(num_SNR, 5);

% Pentru salvarea constelațiilor la praguri importante
constellation_data = struct();

fprintf('\n========================================\n');
fprintf('SIMULARE SNR 0-20 dB (16-QAM)\n');
fprintf('========================================\n');
fprintf('Număr simboluri: %d\n', N);
fprintf('Fading strength: %.1f%%\n', fadingStrength*100);
fprintf('Phase offset: %.2f rad\n', phaseOffset);
fprintf('Prag BER acceptabil: %.1e\n', BER_threshold);
fprintf('Prag BER foarte bun: %.1e\n', BER_threshold_excellent);
fprintf('Prag EVM: %.1f%%\n', EVM_threshold);
fprintf('Prag SNR: %.1f dB\n', SNR_threshold);
fprintf('========================================\n\n');

%% GENERARE SEMNAL TX
rng(42);

data = randi([0 1], N*k, 1);

txSig = qammod(data, M, ...
    'InputType', 'bit', ...
    'UnitAveragePower', true);

% Oversampling și pulse shaping
txUp = upsample(txSig, sps);

rcFilt = rcosdesign(rolloff, span, sps, 'sqrt');

txWave = conv(txUp, rcFilt, 'same');

% Parametri IIR
Fs_oversamp = Fs * sps;

Fc_iir = Fs/2 * (1 + rolloff) * 0.9;

Wn_iir = Fc_iir / (Fs_oversamp/2);

[b_iir, a_iir] = butter(2, Wn_iir);

%% LOOP PRIN NIVELURILE SNR
for idx_snr = 1:num_SNR

    SNR_dB = SNR_range(idx_snr);

    fprintf('Procesare SNR = %d dB... [%d/%d]\n', ...
        SNR_dB, idx_snr, num_SNR);

    %% GENERARE CANAL

    % Fading Rayleigh
    fadingCoeff_wave = ...
        (randn(size(txWave)) + 1j*randn(size(txWave))) / sqrt(2);

    fadingCoeff_wave = ...
        (1-fadingStrength) + fadingStrength*fadingCoeff_wave;

    % Canal pentru semnalul cu oversampling
    rxWave = txWave .* fadingCoeff_wave * exp(1j*phaseOffset);

    rxWave = awgn(rxWave, SNR_dB, 'measured');

    % Canal pentru semnalul la nivel de simboluri
    rxAWGN = awgn(txSig, SNR_dB, 'measured');

    carrierOffset = exp(1j * phaseOffset);

    rxSig = rxAWGN * carrierOffset;

    %% FILTRARE 1: FIR

    rxFilt = conv(rxWave, rcFilt, 'same');

    delay_FIR = span*sps/2;

    rxSampled_FIR = rxFilt(delay_FIR+1:sps:end);

    if length(rxSampled_FIR) < N

        rxFIR = [rxSampled_FIR; ...
            zeros(N-length(rxSampled_FIR), 1)];

    else

        rxFIR = rxSampled_FIR(1:N);

    end

    %% FILTRARE 2: IIR

    rxMatched_IIR = conv(rxWave, rcFilt, 'same');

    rxIIR_wave = filtfilt(b_iir, a_iir, rxMatched_IIR);

    delay_IIR = span*sps/2;

    rxSampled_IIR = rxIIR_wave(delay_IIR+1:sps:end);

    if length(rxSampled_IIR) < N

        rxIIR = [rxSampled_IIR; ...
            zeros(N-length(rxSampled_IIR), 1)];

    else

        rxIIR = rxSampled_IIR(1:N);

    end

    %% FILTRARE 3: LMS

    rxMatched_LMS = conv(rxWave, rcFilt, 'same');

    delay_LMS = span*sps/2;

    rxDownsampled_LMS = rxMatched_LMS(delay_LMS+1:sps:end);

    Nsym_LMS = min(length(rxDownsampled_LMS), N);

    rxDownsampled_LMS = rxDownsampled_LMS(1:Nsym_LMS);

    txAligned_LMS = txSig(1:Nsym_LMS);

    % Compensarea fazei inițiale
    phaseEst_LMS_init = angle(mean( ...
        (rxDownsampled_LMS.^4) .* conj(txAligned_LMS.^4))) / 4;

    rxDownsampled_LMS = ...
        rxDownsampled_LMS * exp(-1j*phaseEst_LMS_init);

    % Normalizare
    rxDownsampled_LMS = ...
        rxDownsampled_LMS / rms(rxDownsampled_LMS);

    % Parametri LMS
    mu = 0.01;
    N_lms = 11;
    N_training = 500;

    w = zeros(N_lms, 1);

    w(ceil(N_lms/2)) = 1;

    y_lms = zeros(Nsym_LMS, 1);

    rxPadded = ...
        [zeros(N_lms-1, 1); rxDownsampled_LMS];

    % Faza de antrenare
    for n = 1:N_training

        idx = n + N_lms - 1;

        x_vec = rxPadded(idx:-1:idx-N_lms+1);

        y_lms(n) = w' * x_vec;

        e = txAligned_LMS(n) - y_lms(n);

        w = w + mu * conj(e) * x_vec;

    end

    % Faza decision-directed
    for n = N_training+1:Nsym_LMS

        idx = n + N_lms - 1;

        x_vec = rxPadded(idx:-1:idx-N_lms+1);

        y_lms(n) = w' * x_vec;

        decision = qammod( ...
            qamdemod(y_lms(n), M, ...
            'UnitAveragePower', true), ...
            M, ...
            'UnitAveragePower', true);

        e = decision - y_lms(n);

        w = w + mu * conj(e) * x_vec;

    end

    y_lms = y_lms / rms(y_lms);

    if length(y_lms) < N

        rxLMS = [y_lms; ...
            zeros(N-length(y_lms), 1)];

    else

        rxLMS = y_lms(1:N);

    end

    %% FILTRARE 4: WAVELET

    rxMatched_Wav = conv(rxWave, rcFilt, 'same');

    waveletName = 'sym4';

    level = 2;

    % Componenta reală
    [C_real, L_real] = ...
        wavedec(real(rxMatched_Wav), level, waveletName);

    sigma_real = ...
        median(abs(C_real(sum(L_real(1:level))+1:end))) / 0.6745;

    thr_real = ...
        sigma_real * sqrt(2*log(length(rxMatched_Wav))) * 0.3;

    C_real_thr = C_real;

    for kk = 1:level

        idx_start = sum(L_real(1:kk)) + 1;

        idx_end = sum(L_real(1:kk+1));

        C_real_thr(idx_start:idx_end) = ...
            wthresh(C_real(idx_start:idx_end), 's', thr_real);

    end

    rxWav_real = ...
        waverec(C_real_thr, L_real, waveletName);

    % Componenta imaginară
    [C_imag, L_imag] = ...
        wavedec(imag(rxMatched_Wav), level, waveletName);

    sigma_imag = ...
        median(abs(C_imag(sum(L_imag(1:level))+1:end))) / 0.6745;

    thr_imag = ...
        sigma_imag * sqrt(2*log(length(rxMatched_Wav))) * 0.3;

    C_imag_thr = C_imag;

    for kk = 1:level

        idx_start = sum(L_imag(1:kk)) + 1;

        idx_end = sum(L_imag(1:kk+1));

        C_imag_thr(idx_start:idx_end) = ...
            wthresh(C_imag(idx_start:idx_end), 's', thr_imag);

    end

    rxWav_imag = ...
        waverec(C_imag_thr, L_imag, waveletName);

    rxWav_denoised = rxWav_real + 1j*rxWav_imag;

    delay_Wav = span*sps/2;

    rxSampled_Wav = rxWav_denoised(delay_Wav+1:sps:end);

    if length(rxSampled_Wav) < N

        rxWavelet = [rxSampled_Wav; ...
            zeros(N-length(rxSampled_Wav), 1)];

    else

        rxWavelet = rxSampled_Wav(1:N);

    end

    %% CALCUL METRICI

    maxShift = 10;

    [EVM_RX, BER_RX, SNR_RX] = ...
        metricsQAM(txSig, rxSig, data, M, k, maxShift);

    [EVM_FIR, BER_FIR, SNR_FIR] = ...
        metricsQAM(txSig, rxFIR, data, M, k, maxShift);

    [EVM_IIR, BER_IIR, SNR_IIR] = ...
        metricsQAM(txSig, rxIIR, data, M, k, maxShift);

    [EVM_LMS, BER_LMS, SNR_LMS] = ...
        metricsQAM(txSig, rxLMS, data, M, k, maxShift);

    [EVM_Wavelet, BER_Wavelet, SNR_Wavelet] = ...
        metricsQAM(txSig, rxWavelet, data, M, k, maxShift);

    %% SALVARE REZULTATE

    EVM_results(idx_snr, :) = ...
        [EVM_RX, EVM_FIR, EVM_IIR, EVM_LMS, EVM_Wavelet];

    BER_results(idx_snr, :) = ...
        [BER_RX, BER_FIR, BER_IIR, BER_LMS, BER_Wavelet];

    SNR_results(idx_snr, :) = ...
        [SNR_RX, SNR_FIR, SNR_IIR, SNR_LMS, SNR_Wavelet];

    %% SALVARE CONSTELAȚII LA PRAGUL BER

    methods = {'RX', 'FIR', 'IIR', 'LMS', 'Wavelet'};

    sigs = {rxSig, rxFIR, rxIIR, rxLMS, rxWavelet};

    for m = 1:5

        field_name = methods{m};

        if ~isfield(constellation_data, field_name) && ...
                BER_results(idx_snr, m) <= BER_threshold

            constellation_data.(field_name).SNR = SNR_dB;

            constellation_data.(field_name).signal = sigs{m};

            constellation_data.(field_name).txSig = txSig;

            constellation_data.(field_name).EVM = ...
                EVM_results(idx_snr, m);

            constellation_data.(field_name).BER = ...
                BER_results(idx_snr, m);

        end

    end

end

fprintf('\n========================================\n');
fprintf('SIMULARE COMPLETĂ!\n');
fprintf('========================================\n\n');

%% TABEL CU REZULTATE

fig_table = uifigure( ...
    'Name', 'Tabel Rezultate SNR vs Performanță', ...
    'Position', [50 50 1500 820]);

method_names = {'RX', 'FIR', 'IIR', 'LMS', 'Wavelet'};

nS = numel(SNR_range);

nM = numel(method_names);

nCols = 1 + 3*nM;

Data = cell(nS+2, nCols);

% Primul rând al antetului
Data(1,1) = {'SNR(dB)'};

c = 2;

for m = 1:nM

    Data(1,c) = method_names(m);

    Data(1,c+1) = {''};

    Data(1,c+2) = {''};

    c = c + 3;

end

% Al doilea rând al antetului
Data(2,1) = {''};

c = 2;

for m = 1:nM

    Data(2,c) = {'BER'};

    Data(2,c+1) = {'EVM(%)'};

    Data(2,c+2) = {'SNR(dB)'};

    c = c + 3;

end

% Datele numerice
Data(3:end,1) = num2cell(SNR_range(:));

c = 2;

for m = 1:nM

    Data(3:end,c) = ...
        num2cell(BER_results(:,m));

    Data(3:end,c+1) = ...
        num2cell(EVM_results(:,m));

    Data(3:end,c+2) = ...
        num2cell(SNR_results(:,m));

    c = c + 3;

end

uit = uitable( ...
    fig_table, ...
    'Data', Data, ...
    'ColumnName', repmat({''}, 1, nCols));

uit.Position = [20 80 1460 720];

uit.ColumnWidth = 'auto';

% Stil pentru antet
s_hdr = uistyle( ...
    'BackgroundColor', [0.9 0.9 0.9], ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

addStyle(uit, s_hdr, 'row', 1);

addStyle(uit, s_hdr, 'row', 2);

% Notă BER
BER_min_res = 1/(N*k);

uilabel(fig_table, ...
    'Text', sprintf(['Notă: BER=0 înseamnă 0 erori observate. ' ...
    'Pentru reprezentare se poate considera BER < ' ...
    '1/(N·k)=%.2e. Praguri: BER acceptabil %.1e, ' ...
    'BER foarte bun %.1e, EVM %.1f%%, SNR %.1f dB.'], ...
    BER_min_res, ...
    BER_threshold, ...
    BER_threshold_excellent, ...
    EVM_threshold, ...
    SNR_threshold), ...
    'Position', [20 40 1460 30]);

% Colorarea tabelului
s_green = uistyle( ...
    'BackgroundColor', [0.6 1 0.6]);

s_red = uistyle( ...
    'BackgroundColor', [1 0.6 0.6]);

for r = 1:nS

    rowUI = r + 2;

    for m = 1:nM

        col_ber = 2 + (m-1)*3;

        col_evm = col_ber + 1;

        col_snr = col_ber + 2;

        if BER_results(r,m) <= BER_threshold

            addStyle(uit, s_green, ...
                'cell', [rowUI, col_ber]);

        end

        if EVM_results(r,m) <= EVM_threshold

            addStyle(uit, s_green, ...
                'cell', [rowUI, col_evm]);

        else

            addStyle(uit, s_red, ...
                'cell', [rowUI, col_evm]);

        end

        if SNR_results(r,m) >= SNR_threshold

            addStyle(uit, s_green, ...
                'cell', [rowUI, col_snr]);

        end

    end

end

%% SETĂRI COMUNE PENTRU CELE TREI GRAFICE

% Grosimea liniilor curbelor
lineWidthCurve = 2.4;

% Grosimea liniilor care marchează pragurile
lineWidthThreshold = 1.8;

% Dimensiunea numerelor de pe axe
fontSizeAxes = 14;

% Dimensiunea denumirilor axelor
fontSizeLabels = 16;

% Dimensiunea titlurilor
fontSizeTitle = 17;

% Dimensiunea textului din legendă
fontSizeLegend = 14;

% Grosimea axelor
axesLineWidth = 1.3;

% Dimensiunea ferestrelor figurilor
figurePosition = [100 100 1250 720];

%% GRAFIC BER vs SNR

figure( ...
    'Name', 'BER vs SNR', ...
    'Position', figurePosition, ...
    'Color', 'w');

% Valorile BER egale cu zero nu pot fi reprezentate direct
% pe scară logaritmică. Pentru afișare, acestea sunt înlocuite
% numai în grafic cu limita minimă 1/(N*k).
BER_min_plot = 1/(N*k);

BER_plot = BER_results;

BER_plot(BER_plot == 0) = BER_min_plot;

semilogy( ...
    SNR_range, ...
    BER_plot(:,1), ...
    '-', ...
    'LineWidth', lineWidthCurve);

hold on;

semilogy( ...
    SNR_range, ...
    BER_plot(:,2), ...
    '--', ...
    'LineWidth', lineWidthCurve);

semilogy( ...
    SNR_range, ...
    BER_plot(:,3), ...
    '-.', ...
    'LineWidth', lineWidthCurve);

semilogy( ...
    SNR_range, ...
    BER_plot(:,4), ...
    ':', ...
    'LineWidth', lineWidthCurve);

semilogy( ...
    SNR_range, ...
    BER_plot(:,5), ...
    '-', ...
    'LineWidth', lineWidthCurve);

yline( ...
    BER_threshold, ...
    '--', ...
    'Prag BER 10^{-4}', ...
    'LineWidth', lineWidthThreshold, ...
    'FontSize', fontSizeAxes, ...
    'LabelHorizontalAlignment', 'right');

yline( ...
    BER_threshold_excellent, ...
    ':', ...
    'Prag BER 10^{-6}', ...
    'LineWidth', lineWidthThreshold, ...
    'FontSize', fontSizeAxes, ...
    'LabelHorizontalAlignment', 'right');

xlabel( ...
    'SNR (dB)', ...
    'FontSize', fontSizeLabels, ...
    'FontWeight', 'bold');

ylabel( ...
    'BER', ...
    'FontSize', fontSizeLabels, ...
    'FontWeight', 'bold');

title( ...
    'Bit Error Rate vs SNR', ...
    'FontSize', fontSizeTitle, ...
    'FontWeight', 'bold');

legend( ...
    'RX', ...
    'FIR', ...
    'IIR', ...
    'LMS', ...
    'Wavelet', ...
    'Prag BER 10^{-4}', ...
    'Prag BER 10^{-6}', ...
    'Location', 'best', ...
    'FontSize', fontSizeLegend);

grid on;

grid minor;

ylim([BER_min_plot/2, 1]);

xlim([min(SNR_range), max(SNR_range)]);

xticks(SNR_range(1):2:SNR_range(end));

ax = gca;

ax.FontSize = fontSizeAxes;

ax.FontWeight = 'bold';

ax.LineWidth = axesLineWidth;

ax.TickDir = 'out';

ax.Box = 'on';

ax.Layer = 'top';

hold off;

%% GRAFIC EVM vs SNR

figure( ...
    'Name', 'EVM vs SNR', ...
    'Position', figurePosition, ...
    'Color', 'w');

plot( ...
    SNR_range, ...
    EVM_results(:,1), ...
    '-', ...
    'LineWidth', lineWidthCurve);

hold on;

plot( ...
    SNR_range, ...
    EVM_results(:,2), ...
    '--', ...
    'LineWidth', lineWidthCurve);

plot( ...
    SNR_range, ...
    EVM_results(:,3), ...
    '-.', ...
    'LineWidth', lineWidthCurve);

plot( ...
    SNR_range, ...
    EVM_results(:,4), ...
    ':', ...
    'LineWidth', lineWidthCurve);

plot( ...
    SNR_range, ...
    EVM_results(:,5), ...
    '-', ...
    'LineWidth', lineWidthCurve);

yline( ...
    EVM_threshold, ...
    '--', ...
    'Prag EVM 12.5%', ...
    'LineWidth', lineWidthThreshold, ...
    'FontSize', fontSizeAxes, ...
    'LabelHorizontalAlignment', 'right');

xlabel( ...
    'SNR (dB)', ...
    'FontSize', fontSizeLabels, ...
    'FontWeight', 'bold');

ylabel( ...
    'EVM (%)', ...
    'FontSize', fontSizeLabels, ...
    'FontWeight', 'bold');

title( ...
    'Error Vector Magnitude vs SNR', ...
    'FontSize', fontSizeTitle, ...
    'FontWeight', 'bold');

legend( ...
    'RX', ...
    'FIR', ...
    'IIR', ...
    'LMS', ...
    'Wavelet', ...
    'Prag EVM', ...
    'Location', 'best', ...
    'FontSize', fontSizeLegend);

grid on;

grid minor;

ylim([0, max(EVM_results(:))*1.1]);

xlim([min(SNR_range), max(SNR_range)]);

xticks(SNR_range(1):2:SNR_range(end));

ax = gca;

ax.FontSize = fontSizeAxes;

ax.FontWeight = 'bold';

ax.LineWidth = axesLineWidth;

ax.TickDir = 'out';

ax.Box = 'on';

ax.Layer = 'top';

hold off;

%% GRAFIC SNR CALCULAT vs SNR SIMULAT

figure( ...
    'Name', 'SNR calculat vs SNR simulat', ...
    'Position', figurePosition, ...
    'Color', 'w');

plot( ...
    SNR_range, ...
    SNR_results(:,1), ...
    '-', ...
    'LineWidth', lineWidthCurve);

hold on;

plot( ...
    SNR_range, ...
    SNR_results(:,2), ...
    '--', ...
    'LineWidth', lineWidthCurve);

plot( ...
    SNR_range, ...
    SNR_results(:,3), ...
    '-.', ...
    'LineWidth', lineWidthCurve);

plot( ...
    SNR_range, ...
    SNR_results(:,4), ...
    ':', ...
    'LineWidth', lineWidthCurve);

plot( ...
    SNR_range, ...
    SNR_results(:,5), ...
    '-', ...
    'LineWidth', lineWidthCurve);

yline( ...
    SNR_threshold, ...
    '--', ...
    'Prag SNR 15 dB', ...
    'LineWidth', lineWidthThreshold, ...
    'FontSize', fontSizeAxes, ...
    'LabelHorizontalAlignment', 'right');

xlabel( ...
    'SNR simulat (dB)', ...
    'FontSize', fontSizeLabels, ...
    'FontWeight', 'bold');

ylabel( ...
    'SNR calculat după procesare (dB)', ...
    'FontSize', fontSizeLabels, ...
    'FontWeight', 'bold');

title( ...
    'Raport Semnal-Zgomot vs SNR simulat', ...
    'FontSize', fontSizeTitle, ...
    'FontWeight', 'bold');

legend( ...
    'RX', ...
    'FIR', ...
    'IIR', ...
    'LMS', ...
    'Wavelet', ...
    'Prag SNR', ...
    'Location', 'best', ...
    'FontSize', fontSizeLegend);

grid on;

grid minor;

xlim([min(SNR_range), max(SNR_range)]);

ylim([min(SNR_results(:))-2, ...
    max(SNR_results(:))*1.1]);

xticks(SNR_range(1):2:SNR_range(end));

ax = gca;

ax.FontSize = fontSizeAxes;

ax.FontWeight = 'bold';

ax.LineWidth = axesLineWidth;

ax.TickDir = 'out';

ax.Box = 'on';

ax.Layer = 'top';

hold off;

%% CONSTELAȚII LA PRAGURI BER IMPORTANTE

if ~isempty(fieldnames(constellation_data))

    fig_const = figure( ...
        'Name', 'Constelații la prag BER', ...
        'Position', [200 200 1200 800]);

    methods_found = fieldnames(constellation_data);

    num_methods = length(methods_found);

    for m = 1:num_methods

        method = methods_found{m};

        data_const = constellation_data.(method);

        subplot(2, 3, m);

        plot( ...
            real(data_const.txSig), ...
            imag(data_const.txSig), ...
            'b.', ...
            'MarkerSize', 3, ...
            'DisplayName', 'TX');

        hold on;

        plot( ...
            real(data_const.signal), ...
            imag(data_const.signal), ...
            'r.', ...
            'MarkerSize', 3, ...
            'DisplayName', method);

        title( ...
            sprintf('%s @ SNR=%d dB\nBER=%.2e, EVM=%.2f%%', ...
            method, ...
            data_const.SNR, ...
            data_const.BER, ...
            data_const.EVM), ...
            'FontSize', 11, ...
            'FontWeight', 'bold');

        xlabel('I (Real)');

        ylabel('Q (Imaginar)');

        legend('Location', 'best');

        grid on;

        axis equal;

        axis([-2 2 -2 2]);

        hold off;

    end

else

    fprintf(['ATENȚIE: Nicio metodă nu a atins pragul ' ...
        'BER <= %.1e în intervalul SNR testat!\n'], ...
        BER_threshold);

end

%% AFIȘARE PRAGURI IMPORTANTE

fprintf('\n========== PRAGURI IMPORTANTE ==========\n');

fprintf('Prag BER acceptabil: %.2e\n', ...
    BER_threshold);

fprintf('Prag BER foarte bun: %.2e\n', ...
    BER_threshold_excellent);

fprintf('Prag EVM pentru 16-QAM: %.1f%%\n', ...
    EVM_threshold);

fprintf('Prag SNR de referință: %.1f dB\n\n', ...
    SNR_threshold);

for m = 1:5

    method = method_names{m};

    % SNR minim pentru BER acceptabil
    idx_ber = find( ...
        BER_results(:,m) <= BER_threshold, ...
        1, ...
        'first');

    if ~isempty(idx_ber)

        fprintf( ...
            '%s - SNR minim pentru BER <= %.1e: %d dB\n', ...
            method, ...
            BER_threshold, ...
            SNR_range(idx_ber));

    else

        fprintf( ...
            '%s - BER > %.1e la toate nivelurile SNR testate!\n', ...
            method, ...
            BER_threshold);

    end

    % SNR minim pentru BER foarte bun
    idx_ber_excellent = find( ...
        BER_results(:,m) <= BER_threshold_excellent, ...
        1, ...
        'first');

    if ~isempty(idx_ber_excellent)

        fprintf( ...
            '%s - SNR minim pentru BER <= %.1e: %d dB\n', ...
            method, ...
            BER_threshold_excellent, ...
            SNR_range(idx_ber_excellent));

    else

        fprintf( ...
            '%s - BER > %.1e la toate nivelurile SNR testate!\n', ...
            method, ...
            BER_threshold_excellent);

    end

    % SNR minim pentru EVM
    idx_evm = find( ...
        EVM_results(:,m) <= EVM_threshold, ...
        1, ...
        'first');

    if ~isempty(idx_evm)

        fprintf( ...
            '%s - SNR minim pentru EVM <= %.1f%%: %d dB\n', ...
            method, ...
            EVM_threshold, ...
            SNR_range(idx_evm));

    else

        fprintf( ...
            '%s - EVM > %.1f%% la toate nivelurile SNR testate!\n', ...
            method, ...
            EVM_threshold);

    end

    % SNR minim pentru SNR calculat
    idx_snr_calc = find( ...
        SNR_results(:,m) >= SNR_threshold, ...
        1, ...
        'first');

    if ~isempty(idx_snr_calc)

        fprintf( ...
            ['%s - SNR simulat minim pentru SNR calculat ' ...
            '>= %.1f dB: %d dB\n'], ...
            method, ...
            SNR_threshold, ...
            SNR_range(idx_snr_calc));

    else

        fprintf( ...
            ['%s - SNR calculat < %.1f dB la toate ' ...
            'nivelurile SNR testate!\n'], ...
            method, ...
            SNR_threshold);

    end

    fprintf('\n');

end