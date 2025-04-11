function H = FiltroEliminaBanda(N, w1, w2)
    % N: N�mero total de puntos
    % w1, w2: l�mites de la banda a eliminar (w1 < w2)
    H = ones(1, N);
    H(w1+1:w2+1) = 0;            % eliminar banda entre w1 y w2
    H(N-w2+1:N-w1+1) = 0;        % eliminar banda sim�trica
end
