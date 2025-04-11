function x_filtrada = FiltroPorNivel(x, L, modo)
    % x           : Se�al original
    % L           : Nivel de umbral (valor absoluto o porcentaje)
    % modo        : 'amplitud' o 'energia'
    %
    % x_filtrada  : Se�al filtrada

    X = fft(x);             % Transformada
    X_mag = abs(X);         % Magnitud
    N = length(X);

    switch lower(modo)
        case 'amplitud'
            % Eliminar componentes con magnitud menor al umbral
            X(X_mag < L) = 0;

        case 'energia'
            % Calcular energ�a total
            energia_total = sum(X_mag.^2);

            % Ordenar y calcular energ�a acumulada
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
            error('Modo inv�lido. Us� ''amplitud'' o ''energia''.')
    end

    x_filtrada = real(ifft(X));  % Se�al filtrada
end
