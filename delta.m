function [d] = delta(r, dr)

    d = zeros(size(r));
    absr = abs(r);

    % Mask accounting for the discrete delta function in matrix form
    m1 = (absr>=0.5*dr) & (absr<=1.5*dr);
    m2 = (~m1) & (absr <= 0.5*dr);
    
    d1 = (5 - 3*absr.*m1./dr -...
        sqrt(-3*(1 - absr.*m1./dr).^2 + 1))./(6*dr);
    d2 = (1 + sqrt(-3*(r.*m2./dr).^2 + 1))./(3*dr);
    
    d = d1.*m1 + d2.*m2;
end

