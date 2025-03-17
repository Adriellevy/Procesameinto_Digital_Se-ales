function [x,Func] = Funcion_seno_Amortiguada(Longitud,Amplitud,Frecuencia,Decaimiento)

    x=0:Longitud;
    Func = Amplitud .*sin(2*pi*Frecuencia*x).* exp(Decaimiento * x);
end

