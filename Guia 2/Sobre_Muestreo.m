function Func = Sobre_Muestreo(Original,Factor)
%SOBRE_MUESTREO 
%   Original: es la se�al que se querra sobremuestrar
%   Factor: ser� cuantas veces se quiere expandir el tama�o de la se�al
%   Se devuelve: [x,y] 
    if Factor <= 0 || mod(Factor, 1) ~= 0
        error('El factor L debe ser un entero positivo.');
    end

    % Longitud de la se�al original
    N = length(Original);

    % Inicializar la se�al sobremuestreada
    Func = zeros(1, N * Factor);

    % Asignar valores
    Func(1:Factor:end) = Original;

    % Visualizaci�n de resultados
    figure;
    stem(0:N-1, Original, 'filled');
    hold on;
    stem(0:(N * Factor - 1), Func, 'r');
    hold off;
    title('Sobremuestreo de la se�al');
    legend('Original', 'Sobremuestreada');

end

