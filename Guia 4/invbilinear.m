function [num_s, den_s] = invbilinear(num_z, den_z, T)
    % INVBILINEAR - Realiza la transformación bilineal inversa (z -> s)
    %
    % [num_s, den_s] = invbilinear(num_z, den_z, T)
    
    % Crear variable simbólica
    syms s z

    % Crear funciones simbólicas en z
    Hz_sym = poly2sym(num_z, z) / poly2sym(den_z, z);

    % Reemplazar z por su equivalente en función de s
    z_s = (1 + s*T/2)/(1 - s*T/2);
    Hs_sym = simplify(subs(Hz_sym, z, z_s));

    % Extraer numerador y denominador
    [num_s_sym, den_s_sym] = numden(Hs_sym);
    num_s = sym2poly(expand(num_s_sym));
    den_s = sym2poly(expand(den_s_sym));

    % Normalizar
    num_s = num_s / den_s(1);
    den_s = den_s / den_s(1);
end
