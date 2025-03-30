function [x,Func] = Funcion_seno_Amortiguada(Longitud,Resolucion,Amplitud,Frecuencia,Decaimiento)

    x=0:Resolucion:Longitud;
    Func = Amplitud .*sin(2*pi*Frecuencia*x).* exp(Decaimiento * x);
end

