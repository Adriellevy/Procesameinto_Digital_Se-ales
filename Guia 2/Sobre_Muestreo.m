function Func = Sobre_Muestreo(Original,Factor)
%SOBRE_MUESTREO 
%   Original: es la señal que se querra sobremuestrar
%   Factor: será cuantas veces se quiere expandir el tamaño de la señal
%   Se devuelve: [x,y] 
    if Factor <= 0 || mod(Factor, 1) ~= 0
        error('El factor L debe ser un entero positivo.');
    end

    % Longitud de la señal original
    N = length(Original);

    % Inicializar la señal sobremuestreada
    Func = zeros(1, N * Factor);

    % Asignar valores
    Func(1:Factor:end) = Original;

    % Visualización de resultados
    figure;
    stem(0:N-1, Original, 'filled');
    hold on;
    stem(0:(N * Factor - 1), Func, 'r');
    hold off;
    title('Sobremuestreo de la señal');
    legend('Original', 'Sobremuestreada');

end

