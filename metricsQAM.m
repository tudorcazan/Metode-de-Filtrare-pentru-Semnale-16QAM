 function [EVM_percent, BER, SNR_dB] = metricsQAM(txSym, rxSym, dataBitsAll, M, k, maxShift)
% metricsQAM: calculează EVM, BER și SNR pentru un semnal QAM
%
% Formule utilizate:
%
% EVM_RMS = sqrt( sum(|rx_corr - tx|^2) / sum(|tx|^2) )
% EVM_%   = 100 * EVM_RMS
%
% BER = numar_biti_eronati / numar_total_biti
%
% SNR = Ps / Pn
% Ps  = mean(|tx|^2)
% Pn  = mean(|rx_corr - tx|^2)
% SNR_dB = 10*log10(Ps/Pn)

    Ntx = length(txSym);
    Nrx = length(rxSym);

    bestEVM = inf;
    bestShift = 0;
    bestStartTx = 1;
    bestStartRx = 1;
    bestNs = min(Ntx, Nrx);
    bestRxCorr = [];
    bestTxAl = [];

    % ==========================================================
    % 1. Căutare aliniere optimă pe simboluri
    %    Alegerea se face după EVM minim, nu după BER.
    % ==========================================================
    for s = -maxShift:maxShift

        if s >= 0
            startTx = 1;
            startRx = 1 + s;
        else
            startTx = 1 - s;
            startRx = 1;
        end

        Ns = min(Ntx - startTx + 1, Nrx - startRx + 1);

        if Ns < 50
            continue;
        end

        tx_al = txSym(startTx:startTx+Ns-1);
        rx_al = rxSym(startRx:startRx+Ns-1);

        % ======================================================
        % 2. Compensare complexă fază + amplitudine
        %    Estimare LS:
        %    alpha = argmin || tx - alpha*rx ||^2
        % ======================================================
        alpha = (rx_al' * tx_al) / (rx_al' * rx_al);
        rx_corr = alpha * rx_al;

        % ======================================================
        % 3. Calcul EVM pentru shift-ul curent
        % ======================================================
        err = rx_corr - tx_al;

        EVM_rms = sqrt(sum(abs(err).^2) / sum(abs(tx_al).^2));

        if EVM_rms < bestEVM
            bestEVM = EVM_rms;
            bestShift = s;
            bestStartTx = startTx;
            bestStartRx = startRx;
            bestNs = Ns;
            bestRxCorr = rx_corr;
            bestTxAl = tx_al;
        end
    end

    % ==========================================================
    % 4. Dacă nu s-a găsit niciun shift valid, folosim aliniere directă
    % ==========================================================
    if isempty(bestRxCorr)

        bestNs = min(Ntx, Nrx);

        bestTxAl = txSym(1:bestNs);
        rx_al = rxSym(1:bestNs);

        alpha = (rx_al' * bestTxAl) / (rx_al' * rx_al);
        bestRxCorr = alpha * rx_al;

        bestStartTx = 1;
    end

    % ==========================================================
    % 5. Calcul final EVM
    % ==========================================================
    err = bestRxCorr - bestTxAl;

    Ps = mean(abs(bestTxAl).^2);
    Pn = mean(abs(err).^2);

    EVM_rms = sqrt(sum(abs(err).^2) / sum(abs(bestTxAl).^2));
    EVM_percent = 100 * EVM_rms;

    % ==========================================================
    % 6. Calcul final SNR
    % ==========================================================
    if Pn > 0
        SNR_dB = 10 * log10(Ps / Pn);
    else
        SNR_dB = inf;
    end

    % ==========================================================
    % 7. Calcul final BER
    % ==========================================================
    startBit = (bestStartTx - 1) * k + 1;
    endBit = startBit + bestNs * k - 1;

    if startBit > length(dataBitsAll)
        BER = NaN;
        return;
    end

    endBit = min(endBit, length(dataBitsAll));
    txBits = dataBitsAll(startBit:endBit);

    rxBits = qamdemod(bestRxCorr, M, ...
        'OutputType', 'bit', ...
        'UnitAveragePower', true);

    Lb = min(length(txBits), length(rxBits));

    if Lb < 1
        BER = NaN;
    else
        numErr = sum(txBits(1:Lb) ~= rxBits(1:Lb));
        BER = numErr / Lb;
    end
end


