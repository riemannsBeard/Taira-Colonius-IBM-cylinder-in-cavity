function [ib] = IBMcylinder(R, nTheta)

    ib.theta = linspace(0, 2*pi*(1 - 1/nTheta), nTheta)';
    ib.r = ib.theta*0 + R;
    
    [ib.xi, ib.eta] = pol2cart(ib.theta, ib.r);
                
    ib.dxi = gradient(ib.xi);
    ib.deta = gradient(ib.eta);
    
end

