function [z, p, k] = my_tf2zp(b, a)
    z = roots(b);
    p = roots(a);
    k = b(1)/a(1);
end