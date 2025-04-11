function x_filtrada = FiltroPorNivel(x, modo, L)
    % x: Señal de entrada (vector de tiempo)
    % modo: 'amplitud' o 'energia'
    % L: umbral (si es 'amplitud', L es el valor de corte; si es 'energia', L es el porcentaje [0-1])
    
    N = length(x);
    X = fft(x);               % Transformada de Fourier
    A = abs(X);               % Magnitud espectral
    
    switch lower(modo)
        case 'amplitud'
            % Umbral por amplitud
            X_filtrado = X .* (A >= L);
            
        case 'energia'
            % Umbral por energía
            [~, idx] = sort(A, 'descend'); % ordenar componentes por magnitud
            energia_total = sum(A.^2);
            energia_acumulada = 0;
            mascara = zeros(size(X));
            
            for i = 1:N
                energia_acumulada = energia_acumulada + A(idx(i))^2;
                mascara(idx(i)) = 1;
                if energia_acumulada / energia_total >= L
                    break;
                end
            end
            
            X_filtrado = X .* mascara;
        otherwise
            error('Modo no válido. Usar "amplitud" o "energia".');
    end
    
    x_filtrada = real(ifft(X_filtrado)); % Señal resultante en el dominio temporal
end
