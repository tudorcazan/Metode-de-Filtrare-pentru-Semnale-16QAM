clc
clear
close all

%%  PARAMETRI SISTEM 
M = 16;                  % Modulație 16-QAM (4 biți/simbol)
N = 250000;              % Număr simboluri pentru statistici robuste
k = log2(M);             % Biți per simbol
Fs = 10e3;               % Rată eșantionare [Hz]
SNR_dB = 15;             % Raport semnal-zgomot [dB]

% Parametri canal    
fadingStrength = 0.1;    % Intensitate fading (0=fără, 1=Rayleigh complet)
phaseOffset    = 0;      % Offset fază [radiani]

% Parametri oversampling și pulse shaping
sps     = 4;             % Samples per symbol (oversampling)
rolloff = 0.25;          % Raised cosine roll-off
span    = 6;             % Filter span în simboluri

% Parametri pentru afișarea constelațiilor
% Metricile sunt calculate pe toate simbolurile, dar constelațiile sunt
% reprezentate pe un subset pentru claritate vizuală.
maxConstPoints = 20000;

%%  GENERARE SEMNAL TX 
rng(42); % Seed pentru reproducibilitate

data  = randi([0 1], N*k, 1);
txSig = qammod(data, M, 'InputType','bit','UnitAveragePower',true);

% Punctele ideale ale constelației 16-QAM
txIdealPoints = unique(txSig);

%%  OVERSAMPLING + PULSE SHAPING (pentru FIR)
txUp   = upsample(txSig, sps);
rcFilt = rcosdesign(rolloff, span, sps, 'sqrt');
txWave = conv(txUp, rcFilt, 'same');

%%  SEMNAL RECEPȚIONAT (RX) - pentru FIR (cu oversampling)

% Fading Rayleigh
fadingCoeff_wave = (randn(size(txWave)) + 1j*randn(size(txWave)))/sqrt(2);
fadingCoeff_wave = (1-fadingStrength) + fadingStrength*fadingCoeff_wave;

% Canal pentru semnal cu oversampling
rxWave = txWave .* fadingCoeff_wave * exp(1j*phaseOffset);
rxWave = awgn(rxWave, SNR_dB, 'measured');

%%  SEMNAL RECEPȚIONAT (RX) - pentru celelalte filtre (fără oversampling)

% Zgomot AWGN
rxAWGN = awgn(txSig, SNR_dB, 'measured');

% Offset de fază
carrierOffset = exp(1j * phaseOffset);

% Semnal recepționat final (RX) pentru IIR, LMS, Wavelet
rxSig = rxAWGN * carrierOffset;

%%  FIGURE 1: SEMNALE ÎN TIMP 
figure('Name','Analiza temporală');

subplot(3,1,1);
plot(1:min(500,N), real(txSig(1:min(500,N))),'r','LineWidth',1.2); 
hold on; 
plot(1:min(500,N), imag(txSig(1:min(500,N))),'b','LineWidth',1.2);
title('Semnal TX (ideal)','FontSize',12,'FontWeight','bold');
xlabel('Index simbol'); 
ylabel('Amplitudine');
legend('I (Real)','Q (Imaginar)','Location','best'); 
grid on; 
xlim([1 min(500,N)]);

subplot(3,1,2);
plot(1:min(500,N), real(rxSig(1:min(500,N))),'r','LineWidth',1.2); 
hold on; 
plot(1:min(500,N), imag(rxSig(1:min(500,N))),'b','LineWidth',1.2);
title(sprintf('Semnal RX (SNR=%d dB, Fading=%.0f%%)', SNR_dB, fadingStrength*100), ...
      'FontSize',12,'FontWeight','bold');
xlabel('Index simbol'); 
ylabel('Amplitudine');
legend('I (Real)','Q (Imaginar)','Location','best'); 
grid on; 
xlim([1 min(500,N)]);

subplot(3,1,3);
error_sig = rxSig - txSig;
plot(1:min(500,N), abs(error_sig(1:min(500,N))),'m','LineWidth',1.2);
title('Eroare RX față de TX (magnitudine)','FontSize',12,'FontWeight','bold');
xlabel('Index simbol'); 
ylabel('|Error|');
grid on; 
xlim([1 min(500,N)]);

%%  FIGURE 2: CONSTELAȚII INIȚIALE 

idxConstInit = randperm(N, min(N, maxConstPoints));

figure('Name','Constelații TX și RX', ...
       'Position',[150 150 1300 450]);

tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile;
scatter(real(txIdealPoints), imag(txIdealPoints), ...
        45, 'b', 'filled');
title('Constelație TX (ideală)','FontSize',12,'FontWeight','bold');
xlabel('I (Real)'); 
ylabel('Q (Imaginar)'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);

nexttile;
scatter(real(rxSig(idxConstInit)), imag(rxSig(idxConstInit)), ...
        3, 'r', 'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeAlpha',0.25);
title(sprintf('Constelație RX (SNR=%d dB)',SNR_dB), ...
      'FontSize',12,'FontWeight','bold');
xlabel('I (Real)'); 
ylabel('Q (Imaginar)'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);

nexttile;
scatter(real(txIdealPoints), imag(txIdealPoints), ...
        45, 'b', 'filled', 'DisplayName','TX ideal');
hold on;
scatter(real(rxSig(idxConstInit)), imag(rxSig(idxConstInit)), ...
        3, 'r', 'filled', ...
        'MarkerFaceAlpha',0.18, ...
        'MarkerEdgeAlpha',0.18, ...
        'DisplayName','RX');
title('Suprapunere TX/RX','FontSize',12,'FontWeight','bold');
xlabel('I (Real)'); 
ylabel('Q (Imaginar)'); 
legend('Location','best'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);
hold off;

%%  FILTRARE 1: FIR - MATCHED FILTER CU OVERSAMPLING 

% Aplicare matched filter
rxFilt = conv(rxWave, rcFilt, 'same');

% Down-sampling după filtrare
delay_FIR     = span*sps/2;
rxSampled_FIR = rxFilt(delay_FIR+1:sps:end);

% Asigurare lungime
if length(rxSampled_FIR) < N
    rxFIR = [rxSampled_FIR; zeros(N - length(rxSampled_FIR), 1)];
else
    rxFIR = rxSampled_FIR(1:N);
end

%%  FILTRARE 2: IIR - MATCHED FILTER + IIR

% PASUL 1: Matched filter SRRC
rxMatched_IIR = conv(rxWave, rcFilt, 'same');

% PASUL 2: Filtru IIR Butterworth lowpass
Fs_oversamp = Fs * sps;
Fc_iir      = Fs/2 * (1 + rolloff) * 0.9;
Wn_iir      = Fc_iir / (Fs_oversamp/2);
[b_iir, a_iir] = butter(2, Wn_iir);

rxIIR_wave = filtfilt(b_iir, a_iir, rxMatched_IIR);

% Down-sampling
delay_IIR     = span*sps/2;
rxSampled_IIR = rxIIR_wave(delay_IIR+1:sps:end);

if length(rxSampled_IIR) < N
    rxIIR = [rxSampled_IIR; zeros(N - length(rxSampled_IIR), 1)];
else
    rxIIR = rxSampled_IIR(1:N);
end

%%  FILTRARE 3: LMS ADAPTIV CU OVERSAMPLING

% PASUL 1: Matched filter SRRC
rxMatched_LMS = conv(rxWave, rcFilt, 'same');

% Down-sampling după matched filter
delay_LMS         = span*sps/2;
rxDownsampled_LMS = rxMatched_LMS(delay_LMS+1:sps:end);

% Aliniere inițială cu TX pentru training
Nsym_LMS          = min(length(rxDownsampled_LMS), N);
rxDownsampled_LMS = rxDownsampled_LMS(1:Nsym_LMS);
txAligned_LMS     = txSig(1:Nsym_LMS);

% Compensare fază inițială
phaseEst_LMS_init = angle(mean((rxDownsampled_LMS.^4) .* conj(txAligned_LMS.^4)))/4;
rxDownsampled_LMS = rxDownsampled_LMS * exp(-1j*phaseEst_LMS_init);

% Normalizare inițială
rxDownsampled_LMS = rxDownsampled_LMS / rms(rxDownsampled_LMS);

% PASUL 2: LMS Adaptiv
mu         = 0.01;          
N_lms      = 11;         
N_training = 5000;   

w = zeros(N_lms, 1);
w(ceil(N_lms/2)) = 1;  

y_lms     = zeros(Nsym_LMS, 1);
e_history = zeros(Nsym_LMS, 1);
w_history = zeros(N_lms, Nsym_LMS);

rxPadded = [zeros(N_lms-1, 1); rxDownsampled_LMS];

% Training
for n = 1:N_training
    idx   = n + N_lms - 1;
    x_vec = rxPadded(idx:-1:idx-N_lms+1);

    y_lms(n) = w' * x_vec;
    e        = txAligned_LMS(n) - y_lms(n);

    e_history(n) = abs(e)^2;
    w = w + mu * conj(e) * x_vec;
    w_history(:,n) = w;
end

% Tracking decision-directed
for n = N_training+1:Nsym_LMS
    idx   = n + N_lms - 1;
    x_vec = rxPadded(idx:-1:idx-N_lms+1);

    y_lms(n) = w' * x_vec;

    decision = qammod(qamdemod(y_lms(n), M, 'UnitAveragePower', true), ...
                      M, 'UnitAveragePower', true);

    e = decision - y_lms(n);

    e_history(n) = abs(e)^2;
    w = w + mu * conj(e) * x_vec;
    w_history(:,n) = w;
end

% Normalizare finală
y_lms = y_lms / rms(y_lms);

% Asigurare lungime
if length(y_lms) < N
    rxLMS = [y_lms; zeros(N - length(y_lms), 1)];
else
    rxLMS = y_lms(1:N);
end

%%  FILTRARE 4: WAVELET DENOISING CU OVERSAMPLING

% PASUL 1: Matched filter SRRC
rxMatched_Wav = conv(rxWave, rcFilt, 'same');

% PASUL 2: Wavelet denoising
waveletName = 'sym4';    
level       = 2;               

[C_real, L_real] = wavedec(real(rxMatched_Wav), level, waveletName);
sigma_real       = median(abs(C_real(sum(L_real(1:level))+1:end))) / 0.6745;
thr_real         = sigma_real * sqrt(2*log(length(rxMatched_Wav))) * 0.3;  

C_real_thr = C_real;

for kk = 1:level
    idx_start = sum(L_real(1:kk)) + 1;
    idx_end   = sum(L_real(1:kk+1));
    C_real_thr(idx_start:idx_end) = wthresh(C_real(idx_start:idx_end), 's', thr_real);
end

rxWav_real = waverec(C_real_thr, L_real, waveletName);

[C_imag, L_imag] = wavedec(imag(rxMatched_Wav), level, waveletName);
sigma_imag       = median(abs(C_imag(sum(L_imag(1:level))+1:end))) / 0.6745;
thr_imag         = sigma_imag * sqrt(2*log(length(rxMatched_Wav))) * 0.3;

C_imag_thr = C_imag;

for kk = 1:level
    idx_start = sum(L_imag(1:kk)) + 1;
    idx_end   = sum(L_imag(1:kk+1));
    C_imag_thr(idx_start:idx_end) = wthresh(C_imag(idx_start:idx_end), 's', thr_imag);
end

rxWav_imag = waverec(C_imag_thr, L_imag, waveletName);

rxWav_denoised = rxWav_real + 1j*rxWav_imag;

% PASUL 3: Down-sampling
delay_Wav     = span*sps/2;
rxSampled_Wav = rxWav_denoised(delay_Wav+1:sps:end);

if length(rxSampled_Wav) < N
    rxWavelet = [rxSampled_Wav; zeros(N - length(rxSampled_Wav), 1)];
else
    rxWavelet = rxSampled_Wav(1:N);
end

%%  FIGURE 3: CONSTELAȚII DUPĂ FILTRARE (NESINCRONIZATE - DOAR VIZUAL)

idxConst = randperm(N, min(N, maxConstPoints));

figure('Name','Constelații după filtrare (brut)');

tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% RX nefiltrat
nexttile;
scatter(real(rxSig(idxConst)), imag(rxSig(idxConst)), ...
        3, 'k', 'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeAlpha',0.25);
title('RX (nefiltrat)','FontSize',12,'FontWeight','bold');
xlabel('I'); 
ylabel('Q'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);

% FIR
nexttile;
scatter(real(rxFIR(idxConst)), imag(rxFIR(idxConst)), ...
        3, 'b', 'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeAlpha',0.25);
title('FIR Matched Filter','FontSize',12,'FontWeight','bold');
xlabel('I'); 
ylabel('Q'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);

% IIR
nexttile;
scatter(real(rxIIR(idxConst)), imag(rxIIR(idxConst)), ...
        3, 'r', 'filled', ...
        'MarkerFaceAlpha',0.20, ...
        'MarkerEdgeAlpha',0.20);
title('IIR filtrat','FontSize',12,'FontWeight','bold');
xlabel('I'); 
ylabel('Q'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);

% LMS
nexttile;
scatter(real(rxLMS(idxConst)), imag(rxLMS(idxConst)), ...
        3, 'g', 'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeAlpha',0.25);
title('LMS adaptiv','FontSize',12,'FontWeight','bold');
xlabel('I'); 
ylabel('Q'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);

% Wavelet
nexttile;
scatter(real(rxWavelet(idxConst)), imag(rxWavelet(idxConst)), ...
        3, 'm', 'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeAlpha',0.25);
title('Wavelet denoising','FontSize',12,'FontWeight','bold');
xlabel('I'); 
ylabel('Q'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);

% TX ideal
nexttile;
scatter(real(txIdealPoints), imag(txIdealPoints), ...
        45, 'c', 'filled');
title('TX (ideal)','FontSize',12,'FontWeight','bold');
xlabel('I'); 
ylabel('Q'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);

%%  SINCRONIZARE TX-RX PE SIMBOLURI

rxSig_sync     = alignSymbolsQAM(txSig, rxSig);
rxFIR_sync     = alignSymbolsQAM(txSig, rxFIR);
rxIIR_sync     = alignSymbolsQAM(txSig, rxIIR);
rxLMS_sync     = alignSymbolsQAM(txSig, rxLMS);
rxWavelet_sync = alignSymbolsQAM(txSig, rxWavelet);

% Lungimi
N_RX  = length(rxSig_sync);
N_FIR = length(rxFIR_sync);
N_IIR = length(rxIIR_sync);
N_LMS = length(rxLMS_sync);
N_WAV = length(rxWavelet_sync);

tx_RX  = txSig(1:N_RX);
tx_FIR = txSig(1:N_FIR);
tx_IIR = txSig(1:N_IIR);
tx_LMS = txSig(1:N_LMS);
tx_WAV = txSig(1:N_WAV);

%%  METRICI OPTIMIZATE (EVM, BER, SNR)

maxShift = 10;  % număr maxim de simboli testați la stânga/dreapta

[EVM_RX,      BER_RX,      SNR_RX]      = metricsQAM(txSig, rxSig,     data, M, k, maxShift);
[EVM_FIR,     BER_FIR,     SNR_FIR]     = metricsQAM(txSig, rxFIR,     data, M, k, maxShift);
[EVM_IIR,     BER_IIR,     SNR_IIR]     = metricsQAM(txSig, rxIIR,     data, M, k, maxShift);
[EVM_LMS,     BER_LMS,     SNR_LMS]     = metricsQAM(txSig, rxLMS,     data, M, k, maxShift);
[EVM_Wavelet, BER_Wavelet, SNR_Wavelet] = metricsQAM(txSig, rxWavelet, data, M, k, maxShift);

fprintf('\n========== REZULTATE METRICI (optimizate) ==========\n');
fprintf('Metodă\t\tEVM(%%)\t\tBER\t\tSNR(dB)\n');

fprintf('RX\t\t%.2f\t\t%.6f\t\t%.2f\n', EVM_RX,      BER_RX,      SNR_RX);
fprintf('FIR\t\t%.2f\t\t%.6f\t\t%.2f\n', EVM_FIR,     BER_FIR,     SNR_FIR);
fprintf('IIR\t\t%.2f\t\t%.6f\t\t%.2f\n', EVM_IIR,     BER_IIR,     SNR_IIR);
fprintf('LMS\t\t%.2f\t\t%.6f\t\t%.2f\n', EVM_LMS,     BER_LMS,     SNR_LMS);
fprintf('Wavelet\t\t%.2f\t\t%.6f\t\t%.2f\n', EVM_Wavelet, BER_Wavelet, SNR_Wavelet);

fprintf('=====================================================\n');

%%  FIGURI SEPARATE: COMPARAȚIE EVM, BER, SNR

methods = {'RX','FIR','IIR','LMS','Wavelet'};

EVM_all = [EVM_RX, EVM_FIR, EVM_IIR, EVM_LMS, EVM_Wavelet];
BER_all = [BER_RX, BER_FIR, BER_IIR, BER_LMS, BER_Wavelet];
SNR_all = [SNR_RX, SNR_FIR, SNR_IIR, SNR_LMS, SNR_Wavelet];

%%  FIGURE 4A: EVM

figure('Name','EVM vs Metodă');

bar(EVM_all, 'FaceColor',[0.2 0.6 0.8]);
set(gca,'XTickLabel',methods,'XTickLabelRotation',45);

ylabel('EVM (%)');
title('Error Vector Magnitude','FontWeight','bold');
grid on;

ylim([0 max(EVM_all)*1.2]);

for i = 1:length(EVM_all)
    text(i, EVM_all(i), sprintf('%.2f', EVM_all(i)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',10);
end

set(gcf,'Position',[300 200 650 450]);

%%  FIGURE 4B: BER

figure('Name','BER vs Metodă');

% Valoarea minimă reprezentabilă statistic pentru BER.
% Dacă BER este 0, înseamnă că nu au fost observate erori în numărul de biți simulați.
% Pentru afișarea pe scară logaritmică, se folosește limita 1/(N*k).
BER_min_plot = 1 / (N * k);

BER_plot = BER_all;
BER_plot(BER_plot == 0) = BER_min_plot;

bar(BER_plot, 'FaceColor',[0.8 0.4 0.2]);
set(gca,'YScale','log');
set(gca,'XTickLabel',methods,'XTickLabelRotation',45);

ylabel('BER');
title('Bit Error Rate','FontWeight','bold');
grid on;

ylim([BER_min_plot/2, max(BER_plot)*10]);

for i = 1:length(BER_all)
    text(i, BER_plot(i)*1.15, sprintf('%.6f', BER_plot(i)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',10);
end

set(gcf,'Position',[300 200 650 450]);

%%  FIGURE 4C: SNR

figure('Name','SNR vs Metodă');

bar(SNR_all, 'FaceColor',[0.4 0.8 0.4]);
set(gca,'XTickLabel',methods,'XTickLabelRotation',45);

ylabel('SNR (dB)');
title('Raport Semnal-Zgomot','FontWeight','bold');
grid on;

ylim([min(SNR_all)-1, max(SNR_all)*1.15]);

for i = 1:length(SNR_all)
    text(i, SNR_all(i), sprintf('%.1f', SNR_all(i)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',10);
end

set(gcf,'Position',[300 200 650 450]);

%%  FIGURE: SNR IMPROVEMENT (față de RX) - fără RX în plot

SNR_in = SNR_RX;  % referința: RX nefiltrat

methods_impr = {'FIR','IIR','LMS','Wavelet'};

SNR_impr = [SNR_FIR - SNR_in, ...
            SNR_IIR - SNR_in, ...
            SNR_LMS - SNR_in, ...
            SNR_Wavelet - SNR_in];

figure('Name','SNR Improvement (față de RX, fără RX)');

bar(SNR_impr);
set(gca,'XTickLabel',methods_impr,'XTickLabelRotation',45);
ylabel('SNR improvement (dB)');
title('Îmbunătățire SNR față de RX (nefiltrat)','FontWeight','bold');
grid on;

yline(0,'k--','LineWidth',1.2);

for i = 1:length(SNR_impr)
    text(i, SNR_impr(i), sprintf('%.2f', SNR_impr(i)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',10);
end

ylim([min(SNR_impr)-1, max(SNR_impr)+1]);

%%  FIGURE 5: ANALIZA SPECTRALĂ

figure('Name','Analiza spectrală');

[pxx_tx,  f_tx]  = pwelch(txSig,       [], [], [], Fs, 'centered');
[pxx_rx,  f_rx]  = pwelch(rxSig_sync,  [], [], [], Fs, 'centered');
[pxx_fir, f_fir] = pwelch(rxFIR_sync,  [], [], [], Fs, 'centered');
[pxx_iir, f_iir] = pwelch(rxIIR_sync,  [], [], [], Fs, 'centered');
[~, ~] = pwelch(rxLMS_sync,  [], [], [], Fs, 'centered');

subplot(2,1,1);
plot(f_tx, 10*log10(pxx_tx), 'k', 'LineWidth',2,'DisplayName','TX ideal'); 
hold on;
plot(f_rx, 10*log10(pxx_rx), 'r', 'LineWidth',1.5,'DisplayName','RX sync');
plot(f_fir,10*log10(pxx_fir),'b','LineWidth',1.5,'DisplayName','FIR sync');
plot(f_iir,10*log10(pxx_iir),'g','LineWidth',1.5,'DisplayName','IIR sync');
xlabel('Frecvență (Hz)'); 
ylabel('PSD (dB/Hz)');
title('Densitate Spectrală de Putere','FontSize',12,'FontWeight','bold');
legend('Location','best'); 
grid on;

%% FIGURE 5: ANALIZA SPECTRALĂ + RĂSPUNSURI FILTRE

figure('Name','Analiza spectrală');

[pxx_tx,  f_tx]  = pwelch(txSig,       [], [], [], Fs, 'centered');
[pxx_rx,  f_rx]  = pwelch(rxSig_sync,  [], [], [], Fs, 'centered');
[pxx_fir, f_fir] = pwelch(rxFIR_sync,  [], [], [], Fs, 'centered');
[pxx_iir, f_iir] = pwelch(rxIIR_sync,  [], [], [], Fs, 'centered');
[pxx_lms, f_lms] = pwelch(rxLMS_sync,  [], [], [], Fs, 'centered');

subplot(3,1,1);
plot(f_tx,  10*log10(pxx_tx),  'k', 'LineWidth',2,   'DisplayName','TX ideal'); 
hold on;
plot(f_rx,  10*log10(pxx_rx),  'r', 'LineWidth',1.5, 'DisplayName','RX sync');
plot(f_fir, 10*log10(pxx_fir), 'b', 'LineWidth',1.5, 'DisplayName','FIR sync');
plot(f_iir, 10*log10(pxx_iir), 'g', 'LineWidth',1.5, 'DisplayName','IIR sync');

xlabel('Frecvență (Hz)'); 
ylabel('PSD (dB/Hz)');
title('Densitate spectrală de putere', ...
      'FontSize',12,'FontWeight','bold');
legend('Location','best'); 
grid on;
hold off;

%% Răspunsurile filtrelor FIR/IIR

[H_fir, W_fir] = freqz(rcFilt, 1, 2048, Fs*sps);
[H_iir, W_iir] = freqz(b_iir, a_iir, 2048, Fs*sps);

% Magnitudine
subplot(3,1,2);
plot(W_fir, 20*log10(abs(H_fir) + eps), ...
     'b', 'LineWidth',2, 'DisplayName','FIR / SRRC');
hold on;
plot(W_iir, 20*log10(abs(H_iir) + eps), ...
     'r', 'LineWidth',2, 'DisplayName','IIR / Butterworth');

xlabel('Frecvență (Hz)');
ylabel('Magnitudine (dB)');
title('Răspunsul în magnitudine al filtrelor FIR și IIR', ...
      'FontSize',12,'FontWeight','bold');
legend('Location','best');
grid on;
hold off;

% Fază
subplot(3,1,3);
plot(W_fir, unwrap(angle(H_fir)), ...
     'b', 'LineWidth',2, 'DisplayName','FIR / SRRC');
hold on;
plot(W_iir, unwrap(angle(H_iir)), ...
     'r', 'LineWidth',2, 'DisplayName','IIR / Butterworth');

xlabel('Frecvență (Hz)');
ylabel('Fază (rad)');
title('Răspunsul în fază al filtrelor FIR și IIR', ...
      'FontSize',12,'FontWeight','bold');
legend('Location','best');
grid on;
hold off;

% set(gcf,'Position',[200 100 950 750]);

%%  FIGURI LMS - ANALIZA FUNCȚIONĂRII FILTRULUI ADAPTIV

% Conversia erorii pătratice în dB
mse_dB = 10*log10(e_history + eps);

% Netezire pentru evidențierea tendinței generale de convergență
winLMS = 500;
mse_dB_smooth = movmean(mse_dB, winLMS);

% Se afișează zona relevantă pentru adaptare, nu toate cele 250.000 de simboluri.
% Astfel se observă clar faza de antrenare și intrarea în regim decision-directed.
maxIterLMS = min(Nsym_LMS, max(30000, 4*N_training));
idxLMS     = 1:maxIterLMS;

% Limite robuste pentru axa Oy, astfel încât vârfurile izolate să nu comprime graficul
ylimLow  = prctile(mse_dB(idxLMS), 1)  - 3;
ylimHigh = prctile(mse_dB(idxLMS), 99) + 3;

if ylimLow >= ylimHigh
    ylimLow  = min(mse_dB(idxLMS)) - 3;
    ylimHigh = max(mse_dB(idxLMS)) + 3;
end

%% FIGURA LMS 1 - Curba de învățare

figure('Name','LMS - Curba de invatare', ...
       'Position',[250 150 950 520]);

hold on;

% Marcarea vizuală a fazei de antrenare
patch([1 N_training N_training 1], ...
      [ylimLow ylimLow ylimHigh ylimHigh], ...
      [0.92 0.92 0.92], ...
      'EdgeColor','none', ...
      'FaceAlpha',0.55, ...
      'DisplayName','Fază de antrenare');

% Marcarea vizuală a fazei decision-directed
patch([N_training maxIterLMS maxIterLMS N_training], ...
      [ylimLow ylimLow ylimHigh ylimHigh], ...
      [0.97 0.97 0.97], ...
      'EdgeColor','none', ...
      'FaceAlpha',0.45, ...
      'DisplayName','Regim decision-directed');

plot(idxLMS, mse_dB(idxLMS), ...
     'Color',[0.65 0.65 0.65], ...
     'LineWidth',0.6, ...
     'DisplayName','Eroare instantanee');

plot(idxLMS, mse_dB_smooth(idxLMS), ...
     'b', ...
     'LineWidth',2.0, ...
     'DisplayName','Eroare mediată');

xline(N_training, 'r--', ...
      'LineWidth',1.8, ...
      'DisplayName','Sfârșit antrenare');

text(N_training/2, ylimHigh-0.08*(ylimHigh-ylimLow), ...
     'Antrenare', ...
     'HorizontalAlignment','center', ...
     'FontWeight','bold', ...
     'Color',[0.25 0.25 0.25]);

text(N_training + (maxIterLMS-N_training)/2, ylimHigh-0.08*(ylimHigh-ylimLow), ...
     'Decision-directed', ...
     'HorizontalAlignment','center', ...
     'FontWeight','bold', ...
     'Color',[0.25 0.25 0.25]);

xlabel('Iterație');
ylabel('MSE (dB)');
title('Curba de învățare a algoritmului LMS', ...
      'FontSize',12,'FontWeight','bold');

legend('Location','best');
grid on;
xlim([1 maxIterLMS]);
ylim([ylimLow ylimHigh]);
set(gca,'FontSize',10);
hold off;

%% FIGURA LMS 2 - Evoluția coeficienților principali

coefCentral = ceil(N_lms/2);
coefToPlot  = unique([coefCentral-2, coefCentral-1, coefCentral, coefCentral+1, coefCentral+2]);
coefToPlot  = coefToPlot(coefToPlot >= 1 & coefToPlot <= N_lms);

% Rărire pentru reprezentare clară, fără pierderea tendinței de adaptare
stepCoef = max(1, floor(maxIterLMS/6000));
idxCoef  = 1:stepCoef:maxIterLMS;

figure('Name','LMS - Evolutia coeficientilor principali', ...
       'Position',[300 180 950 520]);

hold on;

for ii = 1:length(coefToPlot)
    c = coefToPlot(ii);
    
    if c == coefCentral
        lineWidth = 2.2;
    else
        lineWidth = 1.3;
    end
    
    plot(idxCoef, abs(w_history(c, idxCoef)), ...
         'LineWidth',lineWidth, ...
         'DisplayName',sprintf('|w_{%d}|', c));
end

xline(N_training, 'r--', ...
      'LineWidth',1.8, ...
      'DisplayName','Sfârșit antrenare');

xlabel('Iterație');
ylabel('Modul coeficient |w_i|');
title('Evoluția coeficienților principali ai filtrului LMS', ...
      'FontSize',12,'FontWeight','bold');

legend('Location','best');
grid on;
xlim([1 maxIterLMS]);
set(gca,'FontSize',10);
hold off;

%%  FIGURE 8: CONSTELAȚII COMPARATIVE (SINCRONIZATE)

figure('Name','Constelații comparative (sync)', ...
       'Position',[100 100 1400 800]);

tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

idx_RX  = randperm(N_RX,  min(N_RX,  maxConstPoints));
idx_FIR = randperm(N_FIR, min(N_FIR, maxConstPoints));
idx_IIR = randperm(N_IIR, min(N_IIR, maxConstPoints));
idx_LMS = randperm(N_LMS, min(N_LMS, maxConstPoints));
idx_WAV = randperm(N_WAV, min(N_WAV, maxConstPoints));

% TX vs RX
nexttile;
scatter(real(txIdealPoints), imag(txIdealPoints), ...
        35, 'b', 'filled', 'DisplayName','TX ideal');
hold on;
scatter(real(rxSig_sync(idx_RX)), imag(rxSig_sync(idx_RX)), ...
        3, 'r', 'filled', ...
        'MarkerFaceAlpha',0.20, ...
        'MarkerEdgeAlpha',0.20, ...
        'DisplayName','RX sync');
title('TX vs RX','FontWeight','bold'); 
xlabel('I'); 
ylabel('Q'); 
legend('Location','best'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);
hold off;

% RX vs FIR
nexttile;
scatter(real(rxSig_sync(idx_RX)), imag(rxSig_sync(idx_RX)), ...
        3, 'r', 'filled', ...
        'MarkerFaceAlpha',0.18, ...
        'MarkerEdgeAlpha',0.18, ...
        'DisplayName','RX');
hold on;
scatter(real(rxFIR_sync(idx_FIR)), imag(rxFIR_sync(idx_FIR)), ...
        3, 'b', 'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeAlpha',0.25, ...
        'DisplayName','FIR');
title('RX vs FIR (sync)','FontWeight','bold'); 
xlabel('I'); 
ylabel('Q'); 
legend('Location','best'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);
hold off;

% RX vs IIR
nexttile;
scatter(real(rxSig_sync(idx_RX)), imag(rxSig_sync(idx_RX)), ...
        3, 'r', 'filled', ...
        'MarkerFaceAlpha',0.18, ...
        'MarkerEdgeAlpha',0.18, ...
        'DisplayName','RX');
hold on;
scatter(real(rxIIR_sync(idx_IIR)), imag(rxIIR_sync(idx_IIR)), ...
        3, 'g', 'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeAlpha',0.25, ...
        'DisplayName','IIR');
title('RX vs IIR (sync)','FontWeight','bold'); 
xlabel('I'); 
ylabel('Q'); 
legend('Location','best'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);
hold off;

% RX vs LMS
nexttile;
scatter(real(rxSig_sync(idx_RX)), imag(rxSig_sync(idx_RX)), ...
        3, 'r', 'filled', ...
        'MarkerFaceAlpha',0.18, ...
        'MarkerEdgeAlpha',0.18, ...
        'DisplayName','RX');
hold on;
scatter(real(rxLMS_sync(idx_LMS)), imag(rxLMS_sync(idx_LMS)), ...
        3, 'y', 'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeAlpha',0.25, ...
        'DisplayName','LMS');
title('RX vs LMS (sync)', 'FontWeight', 'bold');
xlabel('I'); 
ylabel('Q');
legend('Location','best'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);
hold off;

% RX vs Wavelet
nexttile;
scatter(real(rxSig_sync(idx_RX)), imag(rxSig_sync(idx_RX)), ...
        3, 'r', 'filled', ...
        'MarkerFaceAlpha',0.18, ...
        'MarkerEdgeAlpha',0.18, ...
        'DisplayName','RX');
hold on;
scatter(real(rxWavelet_sync(idx_WAV)), imag(rxWavelet_sync(idx_WAV)), ...
        3, 'c', 'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeAlpha',0.25, ...
        'DisplayName','Wavelet');
title('RX vs Wavelet (sync)','FontWeight','bold'); 
xlabel('I'); 
ylabel('Q'); 
legend('Location','best'); 
grid on; 
axis equal; 
axis([-2 2 -2 2]);
set(gca,'FontSize',10);
hold off;

%%  FIGURE 10: INTERCORELAȚIE (SINCRONIZATĂ)

figure('Name','Intercorelație TX vs Semnale Filtrate (sync)');

[xcorr_RX,      lags_RX]      = xcorr(tx_RX,  rxSig_sync,     'normalized');
[xcorr_FIR,     lags_FIR]     = xcorr(tx_FIR, rxFIR_sync,     'normalized');
[xcorr_IIR,     lags_IIR]     = xcorr(tx_IIR, rxIIR_sync,     'normalized');
[xcorr_LMS,     lags_LMS]     = xcorr(tx_LMS, rxLMS_sync,     'normalized');
[xcorr_Wavelet, lags_Wavelet] = xcorr(tx_WAV, rxWavelet_sync, 'normalized');

subplot(2,3,1);
plot(lags_RX, abs(xcorr_RX),'r','LineWidth',1.2);
title('TX vs RX','FontWeight','bold');
xlabel('Lag'); 
ylabel('|R_{xy}|'); 
grid on; 
xlim([-100 100]);

subplot(2,3,2);
plot(lags_FIR, abs(xcorr_FIR),'b','LineWidth',1.2);
title('TX vs FIR','FontWeight','bold');
xlabel('Lag'); 
ylabel('|R_{xy}|'); 
grid on; 
xlim([-100 100]);

subplot(2,3,3);
plot(lags_IIR, abs(xcorr_IIR),'Color',[0 0.6 0],'LineWidth',1.2);
title('TX vs IIR','FontWeight','bold');
xlabel('Lag'); 
ylabel('|R_{xy}|'); 
grid on; 
xlim([-100 100]);

subplot(2,3,4);
plot(lags_LMS, abs(xcorr_LMS),'m','LineWidth',1.2);
title('TX vs LMS','FontWeight','bold');
xlabel('Lag'); 
ylabel('|R_{xy}|'); 
grid on; 
xlim([-100 100]);

subplot(2,3,5);
plot(lags_Wavelet, abs(xcorr_Wavelet),'c','LineWidth',1.2);
title('TX vs Wavelet','FontWeight','bold');
xlabel('Lag'); 
ylabel('|R_{xy}|'); 
grid on; 
xlim([-100 100]);

subplot(2,3,6);
hold on;
plot(lags_RX, abs(xcorr_RX),'r','LineWidth',1.2,'DisplayName','RX');
plot(lags_FIR, abs(xcorr_FIR),'b','LineWidth',1.2,'DisplayName','FIR');
plot(lags_IIR, abs(xcorr_IIR),'Color',[0 0.6 0],'LineWidth',1.2,'DisplayName','IIR');
plot(lags_LMS, abs(xcorr_LMS),'m','LineWidth',1.2,'DisplayName','LMS');
plot(lags_Wavelet, abs(xcorr_Wavelet),'c','LineWidth',1.2,'DisplayName','Wavelet');
title('Comparație','FontWeight','bold');
xlabel('Lag'); 
ylabel('|R_{xy}|'); 
legend('Location','best'); 
grid on; 
xlim([-100 100]);

%%  CONCLUZIE 

fprintf('\n========== CONCLUZIE ==========\n');

[min_EVM, idx_EVM] = min([EVM_FIR, EVM_IIR, EVM_LMS, EVM_Wavelet]);
[min_BER, idx_BER] = min([BER_FIR, BER_IIR, BER_LMS, BER_Wavelet]);
[max_SNR, idx_SNR] = max([SNR_FIR, SNR_IIR, SNR_LMS, SNR_Wavelet]);

methods_short = {'FIR','IIR','LMS','Wavelet'};

fprintf('Cel mai bun EVM: %s (%.2f%%)\n', methods_short{idx_EVM}, min_EVM);
fprintf('Cel mai bun BER: %s (%.2e)\n', methods_short{idx_BER}, min_BER);
fprintf('Cel mai bun SNR: %s (%.2f dB)\n', methods_short{idx_SNR}, max_SNR);

fprintf('================================\n');

%% Răspunsurile filtrelor FIR/IIR - magnitudine și fază

[H_fir, W_fir] = freqz(rcFilt, 1, 4096, Fs*sps);
[H_iir, W_iir] = freqz(b_iir, a_iir, 4096, Fs*sps);

figure('Name','Răspunsurile filtrelor FIR și IIR', ...
       'Position',[200 100 1050 650]);

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% Magnitudine
nexttile;
plot(W_fir, 20*log10(abs(H_fir) + eps), ...
     'b', 'LineWidth',2, 'DisplayName','FIR / SRRC');
hold on;
plot(W_iir, 20*log10(abs(H_iir) + eps), ...
     'r', 'LineWidth',2, 'DisplayName','IIR / Butterworth');

xlabel('Frecvență (Hz)');
ylabel('Magnitudine (dB)');
title('Răspunsul în magnitudine al filtrelor FIR și IIR', ...
      'FontSize',12,'FontWeight','bold');
legend('Location','best');
grid on;
xlim([0 Fs*sps/2]);
ylim([-150 10]);
set(gca,'FontSize',10);
hold off;

% Fază
nexttile;
plot(W_fir, unwrap(angle(H_fir)), ...
     'b', 'LineWidth',2, 'DisplayName','FIR / SRRC');
hold on;
plot(W_iir, unwrap(angle(H_iir)), ...
     'r', 'LineWidth',2, 'DisplayName','IIR / Butterworth');

xlabel('Frecvență (Hz)');
ylabel('Fază (rad)');
title('Răspunsul în fază al filtrelor FIR și IIR', ...
      'FontSize',12,'FontWeight','bold');
legend('Location','best');
grid on;
xlim([0 Fs*sps/2]);
set(gca,'FontSize',10);
hold off;