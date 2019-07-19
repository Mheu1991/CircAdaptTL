function [AttApproxVeloc,c_re,Z] = CalcApproximateVelocityProfile(c0,Component)
% Calculate lumped parameter elements (i.e. L(alpha), R(alpha)) from
% See Bessems et al. 2007, J Fluid Mech and calculate attenuation constant, wave speed
% and wave impedance as a function of the (characteristic) Womersley number
% M Heusinkveld, Maastricht University, August 2018

global P;

nu        = P.General.Viscosity;              %[m^2/s]
rho       = P.General.rhob;                   %[kg/m3]
omega     = 2*pi* (1/P.General.tCycle);       %[rad/s]

if strcmp(Component,'Tube')
    iComponent = 1:max(size(P.Tube.Name));
    A0         = P.Tube.A0(iComponent);       %[m^2]
    r0         = sqrt(A0/pi);                 %[m]    

elseif strcmp(Component,'ArtVen')
    iComponent = 1:max(size(P.ArtVen.Name));
    A0         = P.ArtVen.A0(:)';             %[m^2]
    r0         = sqrt(A0/pi);                 %[m]
end    

alpha     = sqrt(A0./pi)*sqrt(omega/nu);      %[-] characteristic Womersley number
C         = bsxfun(@rdivide,A0,(rho.*c0.^2)); %[Pa/m^2]

cp = @(alpha) (alpha>sqrt(2)) .* (1+ (sqrt(2)./alpha).*(1-(sqrt(2)./(2*alpha)))) + (alpha<=sqrt(2)).*(3/2);
cq = @(alpha) (alpha>sqrt(2)) .* (alpha./(4*sqrt(2)).*(1-(sqrt(2)./(2*alpha))).^-1) + (alpha<=sqrt(2)).*(1/2);
R0 = @(r0) (8*rho*nu)./(pi*r0.^4);
L0 = @(A0) rho./A0;

%Calculate characteristic Womersley number dependent L and R
R = R0(r0).*(cq(alpha)./(2-cp(alpha)));
L = L0(A0).*(1./(2-cp(alpha)));

AttApproxVeloc = bsxfun(@rdivide,R,(2*sqrt(bsxfun(@rdivide,L,C)))); % R /2*SQRT(L*C)
c_re           = 1./sqrt(bsxfun(@times,L,C));                       % 1/SQRT(L.*C);
Z              = sqrt(bsxfun(@rdivide,L,C));                        % SQRT(L./C);


end