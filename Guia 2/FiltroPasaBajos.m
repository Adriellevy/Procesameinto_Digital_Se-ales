function [Filtro] = FiltroPasaBajos(Longitud,wc)
%FILTROPASABAJOS Summary of this function goes here
%   Detailed explanation goes here
    x = ones(1,Longitud);
    x(wc+1:Longitud-wc)=0;
    figure(2)
    stem(x)
    Filtro = x;
end

