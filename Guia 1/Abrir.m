%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Apertura de señales adquiridas                                          %
%                                                                         %
% function x = Abrir()                                                    %
%                                                                         %
% Dr. Ing. Franco Martin Pessana                                          %
%                                                                         % 
% Facultad de Ingeniería y Ciencias Exactas y Naturales                   %
% Universidad Favaloro                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function x = Abrir()

[Signal,Ubic] = uigetfile('*.txt','Abrir señal discreta');
archivo=[Ubic Signal];
Punt = fopen(archivo,'r');
x = fscanf(Punt,'%f\n');
fclose(Punt);
x = x';