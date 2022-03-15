function xn = normAlong(x, p, dim)
% NORMALONG computes the norm of x along the dimention dim
if dim > ndims(x)
    errorStruct.message = 'dim > ndims(x)';
    error(errorStruct);
end

xn = sqrt(sum(abs(x).^p, dim));
end