function [x, func] = Pulso(Tfinal, Longitud, Amplitud)
    x = 0:Tfinal;

    func = zeros(size(x));

    % Calcular el inicio y el fin del pulso
    inicio = 1;  % Empieza desde el primer índice
    fin = inicio + Longitud;

    % Asegurar que el final no se pase del tiempo total
    fin = min(fin, length(x));

    % Generar el pulso rectangular con la amplitud especificada
    func(inicio:fin) = Amplitud;
    
end