function [xn,xn_retrasada,xm] = Generador_Xm(A,fx,fs,n,M)
xn = A * cos (2*pi*(fx/fs)*n);
xn_retrasada = [zeros(1,M),xn(1:end-M)];
xm = xn .* xn_retrasada;
end