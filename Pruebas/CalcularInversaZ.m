function [ B ] = CalcularInversaZ(x)
    B = ifft(1./fft(x)); 
end

