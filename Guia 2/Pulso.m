function [x, func] = Pulso(Tinicio, Tfinal, InicioPulso, Longitud, Amplitud)
    % Definir el paso (ajustar si es necesario)
    paso = 1;  
    x = Tinicio:paso:Tfinal;  % Vector de tiempo

    func = zeros(size(x)); % Inicializar el vector de funci�n

    % Convertir el inicio del pulso a �ndice
    [~, inicio] = min(abs(x - InicioPulso));  % Encuentra el �ndice m�s cercano a InicioPulso
    [~, fin] = min(abs(x - (InicioPulso + Longitud))); % Encuentra el �ndice para el final del pulso

    % Asegurar que los �ndices sean v�lidos
    inicio = max(inicio, 1);
    fin = min(fin, length(x));

    % Generar el pulso rectangular
    func(inicio:fin) = Amplitud;
end
