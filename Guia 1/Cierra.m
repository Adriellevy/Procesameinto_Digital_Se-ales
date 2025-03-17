function [x, func] = Cierra(L, Amplitud_max, Periodos, Frecuencia)
    % Función para generar una señal diente de sierra
    % L: Número de muestras deseadas
    % A: Amplitud pico de la señal
    % P: Número de períodos en la longitud deseada
    % fT: Frecuencia de muestreo en Hz

    Ts = 1 / Frecuencia;         % Intervalo de muestreo
    t = (0:L-1) * Ts;    % Vector de tiempo discreto
    f0 = Periodos / (L * Ts);   % Frecuencia de la señal

    % Generar la señal diente de sierra
    func = Amplitud_max * sawtooth(2 * pi * f0 * t);
    x = t;
