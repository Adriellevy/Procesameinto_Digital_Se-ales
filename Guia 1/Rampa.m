function [x, func] = Escalon(Longitud,Frecuencia,Pendiente)
%Creacion de una rampa con pendiente variable
%   Detailed explanation goes here
x = 0:1/Frecuencia:Longitud-1;
func = Pendiente*x;

end

