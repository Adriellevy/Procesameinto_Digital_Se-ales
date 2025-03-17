function H = Funcion_Respuesta_discreta(CoefY,CoefX,Longitud)
    h = impz(CoefX,CoefY,Longitud);% Graficar la respuesta al impulso
    stem(0:Longitud-1, h, 'filled')
    xlabel('n')
    ylabel('h[n]')
    title('Respuesta al impulso')
    grid on
end

