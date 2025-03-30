function [x,Func ] = Funcion_Sinc(long,Resolucion,Frecuencia,Desplazamiento )
    x=0:Resolucion:long; 
    %Func = sin(2*pi*Frecuencia .* (x-Desplazamiento))./(2*pi*Frecuencia.*(x-Desplazamiento));
    Func = sinc(2 *pi* Frecuencia .* x-Desplazamiento);
end

