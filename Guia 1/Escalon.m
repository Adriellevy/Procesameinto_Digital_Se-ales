function [x, func] = Escalon(Longitud,Frecuencia)
%Creacion de una funcion escalon
%   Detailed explanation goes here
x = 0:1/Frecuencia:Longitud-1;
func = ones(size(x));

end

