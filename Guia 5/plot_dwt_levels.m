function plot_dwt_levels(x, fs, nivel_max, onda)
% plot_dwt_levels(x, fs, nivel_max, onda)
% ------------------------------------------------------------
% x          ? Se�al de entrada
% fs         ? Frecuencia de muestreo (solo para t�tulo)
% nivel_max  ? Nivel de descomposici�n deseado
% onda       ? Nombre de la ondita madre (ej: 'db4', 'haar', 'sym4', etc.)
%
% Muestra en una sola figura la aproximaci�n y todos los detalles
% correspondientes a la descomposici�n DWT hasta el nivel especificado.

    % Verificar m�ximo nivel permitido
    maxNivel = wmaxlev(length(x), onda);
    if nivel_max > maxNivel
        warning(['Nivel solicitado supera el m�ximo permitido (', num2str(maxNivel), '). Ajustando...']);
        nivel_max = maxNivel;
    end

    % Descomposici�n DWT
    [C, L] = wavedec(x, nivel_max, onda);

    % Aproximaci�n del nivel m�ximo
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
    title(['Aproximaci�n A', num2str(nivel_max)]);
    ylabel(['A', num2str(nivel_max)]);
    grid on;

    for k = 1:nivel_max
        subplot(nivel_max+1, 1, k+1);
        plot(D{k}, 'LineWidth', 1.2);
        title(['Detalle D', num2str(k)]);
        ylabel(['D', num2str(k)]);
        grid on;
    end

    xlabel('�ndice de coeficiente');
    title(['DWT de la se�al Hz � ondita: ', onda]);

end
