function [x, func] = Funcion_senoidal_w(xi,xf, Amplitud_max, FrecAng,Fase)
    % Función para generar una señal cuadrada con ciclo de trabajo variable
    % L: Número de muestras deseadas
    % A: Amplitud pico de la señal

    t = (xi:xf);    % Vector de tiempo discreto

    % Generar la señal cuadrada
    func = Amplitud_max * sin(FrecAng * t + Fase);
    x = t;
end