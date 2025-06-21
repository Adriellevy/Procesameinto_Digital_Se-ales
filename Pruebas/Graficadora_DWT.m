function Graficadora_DWT(x, fs, nivel_max, onda)
% Graficadora_DWT(x, fs, nivel_max, onda)
% ------------------------------------------------------------
% x          ? Señal de entrada
% fs         ? Frecuencia de muestreo
% nivel_max  ? Nivel de descomposición deseado
% onda       ? Nombre de la ondita madre (ej: 'db4', 'haar', 'sym4', etc.)
%
% Muestra la aproximación y los detalles de la DWT, indicando
% las bandas de frecuencia asociadas a cada nivel.

    % Verificar nivel permitido
    maxNivel = wmaxlev(length(x), onda);
    if nivel_max > maxNivel
        warning(['Nivel solicitado supera el máximo permitido (', num2str(maxNivel), '). Ajustando...']);
        nivel_max = maxNivel;
    end

    % Descomposición
    [C, L] = wavedec(x, nivel_max, onda);
    A = appcoef(C, L, onda, nivel_max);
    D = cell(1, nivel_max);
    for k = 1:nivel_max
        D{k} = detcoef(C, L, k);
    end

    % Parámetro para agrupar subplots por figura
    plots_per_figure = 6;
    total_plots = nivel_max + 1; % Aproximación + detalles
    num_figures = ceil(total_plots / plots_per_figure);

    % Graficar en varias figuras si es necesario
    plot_index = 1;
    for fig = 1:num_figures
        figure('Name', ['DWT - ', onda, ' - Fig ', num2str(fig)], 'NumberTitle', 'off');
        for sp = 1:min(plots_per_figure, total_plots - plots_per_figure * (fig - 1))
            subplot(plots_per_figure, 1, sp);

            if plot_index == 1
                % Aproximación
                f_max = fs / 2^nivel_max;
                plot(A, 'LineWidth', 1.2);
                title(sprintf('Aproximación A%d ~ [0 - %.1f] Hz', nivel_max, f_max));
                ylabel(['A', num2str(nivel_max)]);
            else
                % Detalles
                idx = plot_index - 1;
                f_low = fs / 2^(idx+1);
                f_high = fs / 2^idx;
                plot(D{idx}, 'LineWidth', 1.2);
                title(sprintf('Detalle D%d ~ [%.1f - %.1f] Hz', idx, f_low, f_high));
                ylabel(['D', num2str(idx)]);
            end

            grid on;
            plot_index = plot_index + 1;
        end

        % Pausar entre figuras si hay más
        if fig < num_figures
            disp('Presioná una tecla para continuar a la siguiente figura...');
            pause;
        end
    end

    xlabel('Índice de coeficiente');
end
