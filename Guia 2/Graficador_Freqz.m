function Graficador_Freqz(H, w)
    % Graficador_Funciones_Complejas - Dibuja 4 gráficas (modulo, fase, real, imag)
    % a partir de la respuesta en frecuencia H y el vector de frecuencias w.
    %
    % Inputs:
    %   H - vector de respuesta en frecuencia (complejo), salida de freqz
    %   w - vector de frecuencias correspondiente (en radianes por muestra)

    % Verificación básica
    if length(H) ~= length(w)
        error('El tamaño de H y w debe coincidir');
    end

    w = w/pi;
    % Calcular componentes
    modulo = abs(H);
    fase = angle(H);
    parte_real = real(H);
    parte_imag = imag(H);

    % Graficar
    figure;

    subplot(2,2,1);
    plot(w, modulo, 'LineWidth', 1.5);
    grid on;
    title('Módulo |H(w)|');
    xlabel('\omega [rad/muestra]'); ylabel('|H(w)|');

    subplot(2,2,2);
    plot(w, fase, 'LineWidth', 1.5);
    grid on;
    title('Fase arg(H(w))');
    xlabel('\omega [rad/muestra]'); ylabel('Fase [rad]');

    subplot(2,2,3);
    plot(w, parte_real, 'LineWidth', 1.5);
    grid on;
    title('Parte Real Re(H(w))');
    xlabel('\omega [rad/muestra]'); ylabel('Re(H(w))');

    subplot(2,2,4);
    plot(w, parte_imag, 'LineWidth', 1.5);
    grid on;
    title('Parte Imaginaria Im(H(w))');
    xlabel('\omega [rad/muestra]'); ylabel('Im(H(w))');
end
