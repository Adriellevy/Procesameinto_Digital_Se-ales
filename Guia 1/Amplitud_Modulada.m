function [x,Func] = Amplitud_Modulada(Longitud,AmplitudMensajera,ModuloMensajera,frecuenciaMensajera,frecuenciaSenial)

    x = 0:Longitud; 
    Func = AmplitudMensajera*(1+ModuloMensajera*sin(2*pi*frecuenciaMensajera*x(:))).*cos(2*pi*frecuenciaSenial*x(:));

end

