function Graficador_Fft(H_fft, Titulo)
    % Graficador_FFT - Dibuja 4 gráficas (módulo, fase, parte real, imaginaria)
    % a partir de la salida de la transformada FFT, solo en el intervalo [0, pi].
    %
    % Inputs:
    %   H_fft - vector complejo de la FFT de la señal
    %   Titulo - (opcional) título de la figura

    N = length(H_fft);
    M = floor(N/2) + 1; % Solo parte positiva del espectro
    H = H_fft(1:M);     % Tomar hasta pi
    w = linspace(0, pi, M); % Eje de frecuencia en rad/muestra

    w_norm = w / pi;  % Para graficar en [0,1]

    % Calcular componentes
    modulo = abs(H);
    fase = angle(H);
    parte_real = real(H);
    parte_imag = imag(H);

    % Crear figura con o sin título
    if nargin == 2
        figure('Name', Titulo);
    else
        figure;
    end

    subplot(2,2,1);
    plot(w_norm, modulo, 'LineWidth', 1.5);
    grid on;
    title('Módulo |H(\omega)|');
    xlabel('\omega / \pi'); ylabel('|H(\omega)|');

    subplot(2,2,2);
    plot(w_norm, fase, 'LineWidth', 1.5);
    grid on;
    title('Fase arg(H(\omega))');
    xlabel('\omega / \pi'); ylabel('Fase [rad]');

    subplot(2,2,3);
    plot(w_norm, parte_real, 'LineWidth', 1.5);
    grid on;
    title('Parte Real Re(H(\omega))');
    xlabel('\omega / \pi'); ylabel('Re(H(\omega))');

    subplot(2,2,4);
    plot(w_norm, parte_imag, 'LineWidth', 1.5);
    grid on;
    title('Parte Imaginaria Im(H(\omega))');
    xlabel('\omega / \pi'); ylabel('Im(H(\omega))');
end
