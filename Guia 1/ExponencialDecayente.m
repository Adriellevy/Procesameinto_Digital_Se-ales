function [x,Func] = ExponencialDecayente(Longitud,Amplitud,Desplazamiento,Decaimiento)

    x=0:Longitud; 
    % Creamos la parte de ceros hasta el desplazamiento
    ParteCeros = zeros(1, Desplazamiento);
    
    % Generamos la parte exponencial a partir del desplazamiento
    ParteExponencial = Amplitud .* exp(Decaimiento * x(1:end-Desplazamiento));

    % Unimos ambas partes
    Func = [ParteCeros, ParteExponencial];
    %Func=;

end

