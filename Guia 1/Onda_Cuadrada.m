function [x, func] = Onda_Cuadrada(L, Amplitud_max, Periodos,Frecuencia, dutyCycle )
    % Función para generar una señal cuadrada con ciclo de trabajo variable
    % L: Número de muestras deseadas
    % A: Amplitud pico de la señal
    % P: Número de períodos en la longitud deseada
    % dutyCycle: Ciclo de trabajo en porcentaje (0-100)
    % fT: Frecuencia de muestreo en Hz

    Ts = 1 / Frecuencia;         % Intervalo de muestreo
    t = (0:L-1) * Ts;    % Vector de tiempo discreto
    f0 = Periodos / (L * Ts);   % Frecuencia de la señal

    % Generar la señal cuadrada
    func = Amplitud_max * square(2 * pi * f0 * t, dutyCycle);
    x = t;
end