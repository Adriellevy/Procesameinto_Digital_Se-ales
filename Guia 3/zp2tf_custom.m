function [b, a] = zp2tf_custom(z, p, k)
    b = k * poly(z);
    a = poly(p);
end