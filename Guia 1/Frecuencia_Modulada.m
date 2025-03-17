function [x, Func] = Frecuencia_Modulada(A, fc, fm, Am, T, Fs)
    % A: Amplitud de la señal portadora
    % fc: Frecuencia de la portadora (Hz)
    % fm: Frecuencia de la moduladora (Hz)
    % Am: Amplitud de la moduladora
    % T: Duración total de la señal (segundos)
    % Fs: Frecuencia de muestreo (Hz)

    % Definimos el tiempo
    x = 0:1/Fs:T;

    % Calculamos la señal FM
    Func = A * cos(2 * pi * fc * x + Am * sin(2 * pi * fm * x));
end