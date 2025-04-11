function x_filtrada = FiltroPorNivel(x, L, modo)
    % x           : Señal original
    % L           : Nivel de umbral (valor absoluto o porcentaje)
    % modo        : 'amplitud' o 'energia'
    %
    % x_filtrada  : Señal filtrada

    X = fft(x);             % Transformada
    X_mag = abs(X);         % Magnitud
    N = length(X);

    switch lower(modo)
        case 'amplitud'
            % Eliminar componentes con magnitud menor al umbral
            X(X_mag < L) = 0;

        case 'energia'
            % Calcular energía total
            energia_total = sum(X_mag.^2);

            % Ordenar y calcular energía acumulada
            energia = X_mag.^2;
            [~, orden] = sort(energia, 'descend');
            energia_acumulada = 0;
            X_nueva = zeros(size(X));

            for i = 1:N
                idx = orden(i);
                energia_acumulada = energia_acumulada + energia(idx);
                X_nueva(idx) = X(idx);
                if energia_acumulada >= (L/100) * energia_total
                    break;
                end
            end

            X = X_nueva;

        otherwise
            error('Modo inválido. Usá ''amplitud'' o ''energia''.')
    end

    x_filtrada = real(ifft(X));  % Señal filtrada
end
