function r = Correlacion(x, y)
    % Si solo se pasa una señal, se hace autocorrelación
    if nargin < 2
        y = x;
    end

    % Asegurar que ambas señales tengan la misma longitud
    N = length(x) + length(y) - 1;
    X = fft(x, N);
    Y = fft(y, N);
    Y = conj(Y);
    % Graficador_Fft(X,'Original')
    % Graficador_Fft(Y,'Conjunta')

    % Correlación cruzada: IFFT del producto de X por conjugado de Y
    R = X .*Y;
    r = ifft(R);

end
