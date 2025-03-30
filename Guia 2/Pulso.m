function [x, func] = Pulso(Tinicio, Tfinal, InicioPulso, Longitud, Amplitud)
    % Definir el paso (ajustar si es necesario)
    paso = 1;  
    x = Tinicio:paso:Tfinal;  % Vector de tiempo

    func = zeros(size(x)); % Inicializar el vector de función

    % Convertir el inicio del pulso a índice
    [~, inicio] = min(abs(x - InicioPulso));  % Encuentra el índice más cercano a InicioPulso
    [~, fin] = min(abs(x - (InicioPulso + Longitud))); % Encuentra el índice para el final del pulso

    % Asegurar que los índices sean válidos
    inicio = max(inicio, 1);
    fin = min(fin, length(x));

    % Generar el pulso rectangular
    func(inicio:fin) = Amplitud;
end
