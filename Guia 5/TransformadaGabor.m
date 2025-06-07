function [GaborTransform, t_new, f_new] = TransformadaGabor(x, fs, Alfa, ModulacionFrecuencial)
    N = length(x);
    Ts = 1 / fs;
    t = 0:Ts:(N - 1) * Ts;
    T0 = N * Ts;
    At = T0 / ModulacionFrecuencial;

    GaborTransform = zeros(ModulacionFrecuencial + 1, N); 

    for k = 1:(ModulacionFrecuencial + 1)
        ventana = exp(-(t - (k - 1) * At).^2 / (4 * Alfa)) / (2 * sqrt(pi * Alfa));
        senial_filtrada = x .* ventana;
        GaborTransform(k, :) = abs(fft(senial_filtrada));
    end

    % Solo el medio espectro
    GaborTransform = GaborTransform(:, 1:floor(N/2));
    t_new = 0:At:T0;
    f_new = linspace(0, fs/2, floor(N/2));
end
