clc
clear
close all

%% FIGURI SEPARATE:
% Reprezentarea în timp, în frecvență și în planul I/Q
% pentru semnalul I/Q 16-QAM generat

%% PARAMETRI SISTEM

M = 16;                  % Modulație 16-QAM
N = 250000;              % Număr simboluri
k = log2(M);             % Biți per simbol
Fs = 10e3;               % Rată de eșantionare simbolică [Hz]

% Parametri oversampling și pulse shaping
sps     = 4;             % Samples per symbol
rolloff = 0.25;          % Factor roll-off
span    = 6;             % Lungimea filtrului în simboluri

% Frecvența efectivă după oversampling
Fs_oversamp = Fs * sps;

% Parametri pentru afișare
nTimePlot = 800;         % Eșantioane afișate în timp
maxConstPoints = 20000;  % Puncte afișate în constelație

%% GENERARE SEMNAL I/Q 16-QAM

rng(42);

data = randi([0 1], N*k, 1);

txSig = qammod(data, M, ...
    'InputType','bit', ...
    'UnitAveragePower',true);

%% OVERSAMPLING + PULSE SHAPING

txUp = upsample(txSig, sps);

rcFilt = rcosdesign(rolloff, span, sps, 'sqrt');

txWave = conv(txUp, rcFilt, 'same');

%% SELECTARE PUNCTE PENTRU CONSTELAȚIE

idxConst = randperm(N, min(N, maxConstPoints));

%% ============================================================
% FIGURA 1: REPREZENTAREA ÎN TIMP A SEMNALULUI I/Q
% ============================================================

figure('Name','Semnal IQ generat - Domeniul timp', ...
       'Position',[250 180 900 500]);

nPlot = min(nTimePlot, length(txWave));

plot(1:nPlot, real(txWave(1:nPlot)), ...
     'b', ...
     'LineWidth',1.2, ...
     'DisplayName','Componenta I');
hold on;

plot(1:nPlot, imag(txWave(1:nPlot)), ...
     'r', ...
     'LineWidth',1.2, ...
     'DisplayName','Componenta Q');

xlabel('Index eșantion');
ylabel('Amplitudine');
title('Reprezentarea în timp a semnalului I/Q generat', ...
      'FontSize',12, ...
      'FontWeight','bold');

legend('Location','best');
grid on;
xlim([1 nPlot]);

set(gca,'FontSize',10);
hold off;



%% FIGURA 2: REPREZENTAREA ÎN FRECVENȚĂ A SEMNALULUI I/Q


% PSD calculată cu fereastră mai mare pentru o reprezentare mai netedă
windowPSD   = hamming(4096);
noverlapPSD = round(0.75 * length(windowPSD));
nfftPSD     = 8192;

[pxx_tx, f_tx] = pwelch(txWave, ...
                        windowPSD, ...
                        noverlapPSD, ...
                        nfftPSD, ...
                        Fs_oversamp, ...
                        'centered');

% Conversie frecvență în kHz
f_kHz = f_tx / 1e3;

% PSD în dB/Hz
PSD_dB = 10*log10(pxx_tx + eps);

% Netezire ușoară pentru aspect mai curat
PSD_dB_smooth = movmean(PSD_dB, 7);

figure('Name','Semnal IQ generat - Domeniul frecvență', ...
       'Position',[300 200 900 500]);

plot(f_kHz, PSD_dB_smooth, ...
     'k', ...
     'LineWidth',1.4);

xlabel('Frecvență (kHz)');
ylabel('PSD (dB/Hz)');
title('Reprezentarea în frecvență a semnalului I/Q', ...
      'FontSize',12, ...
      'FontWeight','bold');

grid on;

% Se afișează zona utilă mai clar, fără scala ×10^4
xlim([-10 10]);

% Limită verticală mai potrivită pentru Word
ylim([max(PSD_dB_smooth)-70, max(PSD_dB_smooth)+5]);

set(gca,'FontSize',10);




%% FIGURA 3: DIAGRAMA I/Q - CONSTELAȚIA 16-QAM


figure('Name','Semnal IQ generat - Diagramă IQ', ...
       'Position',[350 220 650 600]);

scatter(real(txSig(idxConst)), ...
        imag(txSig(idxConst)), ...
        6, ...
        'b', ...
        'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeAlpha',0.25);

xlabel('I');
ylabel('Q');
title('Diagrama I/Q a semnalului 16-QAM generat', ...
      'FontSize',12, ...
      'FontWeight','bold');

grid on;
axis equal;
axis([-1.5 1.5 -1.5 1.5]);

set(gca,'FontSize',10);

