function [x, Func] = Frecuencia_Modulada(A, fc, fm, Am, T, Fs)
    % A: Amplitud de la se�al portadora
    % fc: Frecuencia de la portadora (Hz)
    % fm: Frecuencia de la moduladora (Hz)
    % Am: Amplitud de la moduladora
    % T: Duraci�n total de la se�al (segundos)
    % Fs: Frecuencia de muestreo (Hz)

    % Definimos el tiempo
    x = 0:1/Fs:T;

    % Calculamos la se�al FM
    Func = A * cos(2 * pi * fc * x + Am * sin(2 * pi * fm * x));
end