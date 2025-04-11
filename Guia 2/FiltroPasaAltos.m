function H = FiltroPasaAltos(N, w1)
    % N: Número total de puntos
    % w1: Frecuencia de corte
    H = zeros(1, N);
    H(w1+1:N-w1+1) = 1;  % pasa desde w1 hasta N-w1
end
