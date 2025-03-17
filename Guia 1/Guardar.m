%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% AGuardado de Señales Procesadas                                         %
%                                                                         %
% function x = Guardar(x)                                                 %
%                                                                         %
% Dr. Ing. Franco Martin Pessana                                          %
%                                                                         % 
% Facultad de Ingeniería y Ciencias Exactas y Naturales                   %
% Universidad Favaloro                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function Guardar(y)

[Signal,Ubic]=uiputfile('*.txt','Guardar Señal Resultante');
Signal=[Signal '.txt'];
archivo=[Ubic Signal];
Punt=fopen(archivo,'w');
fprintf(Punt,'%4.2f\n',y);
fclose(Punt);
return


