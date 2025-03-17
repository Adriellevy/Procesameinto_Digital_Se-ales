function H = Funcion_Respuesta_discreta(CoefX,CoefY,Longitud)
    H = impz(CoefX,CoefY,Longitud);% Graficar la respuesta al impulso
    stem(0:Longitud-1, H, 'filled')
    xlabel('n')
    ylabel('h[n]')
    title('Respuesta al impulso')
    grid on
end

