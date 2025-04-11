function H = FiltroPasaBanda(N, w1, w2)
    % N: Número total de puntos
    % w1, w2: límites de la banda que se desea pasar
    H = zeros(1, N);
    H(w1+1:w2+1) = 1;             % banda directa
    H(N-w2+1:N-w1+1) = 1;         % banda simétrica
end
