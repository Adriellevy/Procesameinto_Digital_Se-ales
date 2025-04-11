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

    %Func = interp(Func,Factor);
    %Func = Func(1:Factor:end);
    %Filtro = FiltroPasaBajos(N*Factor,pi/Factor);
    %FuncionEnFourrier = fft(Func);
    %Func = Filtro.*FuncionEnFourrier; 
    %Func = ifft(Func);

end

