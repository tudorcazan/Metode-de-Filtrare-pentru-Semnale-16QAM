# Tehnici de reducere a zgomotului în semnale 16-QAM

Acest repository conține codurile MATLAB utilizate pentru realizarea părții aplicative a proiectului de diplomă „Tehnici de reducere a zgomotului în semnale de radiocomunicații”.

Simularea analizează comparativ performanța unor metode de reducere a zgomotului aplicate unui semnal modulat 16-QAM:

- filtrare FIR;
- filtrare IIR;
- filtrare adaptivă LMS;
- denoising bazat pe transformata wavelet.

Evaluarea performanței este realizată pe baza metricilor Error Vector Magnitude, Bit Error Rate și Signal-to-Noise Ratio.

## Fișierele repository-ului

- `semnal_IQ_generat.m` — generarea și reprezentarea în timp, frecvență și plan I/Q a semnalului 16-QAM;
- `cod_qam.m` — simularea principală, modelarea canalului, aplicarea metodelor FIR, IIR, LMS și Wavelet și generarea rezultatelor comparative;
- `antrenare_lms.m` — analiza influenței numărului de simboluri asupra performanței metodei LMS;
- `diagrama_LMS.m` — generarea schemei de funcționare a filtrului adaptiv LMS;
- `SNR_0_20.m` — analiza performanței sistemului pentru valori ale SNR cuprinse între 0 și 20 dB;
- `alignSymbolsQAM.m` — funcție pentru sincronizarea simbolurilor transmise și recepționate;
- `metricsQAM.m` — funcție pentru calculul metricilor EVM, BER și SNR.

## Cerințe software

Codurile au fost dezvoltate și testate în MATLAB și utilizează funcții din:

- Communications Toolbox;
- Signal Processing Toolbox;
- Wavelet Toolbox.
