function [ Func ] = Sub_Muestreo(Original,Factor)
%SUB_MUESTREO Summary of this function goes here
%   Detailed explanation goes here

    Func = Original(1:Factor:end);
end

