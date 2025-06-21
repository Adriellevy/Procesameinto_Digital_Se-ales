function plot_dwt_levels(x, fs, nivel_max, onda)
% plot_dwt_levels(x, fs, nivel_max, onda)
% ------------------------------------------------------------
% x          ? Señal de entrada
% fs         ? Frecuencia de muestreo (solo para título)
% nivel_max  ? Nivel de descomposición deseado
% onda       ? Nombre de la ondita madre (ej: 'db4', 'haar', 'sym4', etc.)
%
% Muestra en una sola figura la aproximación y todos los detalles
% correspondientes a la descomposición DWT hasta el nivel especificado.

    % Verificar máximo nivel permitido
    maxNivel = wmaxlev(length(x), onda);
    if nivel_max > maxNivel
        warning(['Nivel solicitado supera el máximo permitido (', num2str(maxNivel), '). Ajustando...']);
        nivel_max = maxNivel;
    end

    % Descomposición DWT
    [C, L] = wavedec(x, nivel_max, onda);

    % Aproximación del nivel máximo
    A = appcoef(C, L, onda, nivel_max);

    % Extraer detalles
    D = cell(1, nivel_max);
    for k = 1:nivel_max
        D{k} = detcoef(C, L, k);
    end

    % Graficar
    figure('Name', ['DWT - ', onda, ' - Nivel ', num2str(nivel_max)], 'NumberTitle', 'off');
    
    subplot(nivel_max+1, 1, 1);
    plot(A, 'LineWidth', 1.2);
    title(['Aproximación A', num2str(nivel_max)]);
    ylabel(['A', num2str(nivel_max)]);
    grid on;

    for k = 1:nivel_max
        subplot(nivel_max+1, 1, k+1);
        plot(D{k}, 'LineWidth', 1.2);
        title(['Detalle D', num2str(k)]);
        ylabel(['D', num2str(k)]);
        grid on;
    end

    xlabel('Índice de coeficiente');
    title(['DWT de la señal Hz – ondita: ', onda]);

end
