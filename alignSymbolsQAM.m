function rx_sync = alignSymbolsQAM(txSym, rxSym)

    % Cross-correlation doar între simboluri (mult mai stabil)
    [c, lags] = xcorr(rxSym, txSym);

    [~, idx] = max(abs(c));
    delay = lags(idx);

    % Ajustare delay
    if delay < 0
        rx_sync = rxSym(-delay+1:end);
    else
        rx_sync = [zeros(delay,1); rxSym];
    end

    % Trunchiere finală
    L = min(length(rx_sync), length(txSym));
    rx_sync = rx_sync(1:L);

end
