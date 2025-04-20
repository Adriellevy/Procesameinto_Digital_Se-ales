function [yFiltrada, CoeficientesFiltro] = FiltroPasaBajosFaseCero(x, Orden, wc)
% FiltroPasaBajos aplica un filtro pasa bajos de fase cero usando filtfilt
%   x: señal de entrada
%   Orden: orden del filtro FIR (número de coeficientes - 1)
%   wc: frecuencia de corte normalizada (entre 0 y 1, donde 1 = Nyquist)

    % Diseño del filtro FIR pasa bajos
    CoeficientesFiltro = fir1(Orden, wc,'low');  % Filtro FIR pasa bajos de fase lineal

    % Aplicación del filtro con fase cero usando filtfilt
    yFiltrada = filtfilt(CoeficientesFiltro, 1, x);
end
