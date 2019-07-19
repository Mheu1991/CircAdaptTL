function ArtVenV2p
%function ArtVenV2p
% Art/Ven hemodynamics with peripheral resistance in between
% Flow waves enter Art and Ven.
% volume V-> transmural pressure pTrans(V) and wave impedance Z(V)
% Theo Arts, Maastricht University, Oct 12, 2012

global P;

ArtVen = P.ArtVen; % ArtVen structure
rhob   = P.General.rhob;
% indices referring to related cavities and walls
iCavity= [ArtVen.iCavity; ArtVen.iCavity+1]; iCavity=iCavity(:);

V      = P.Cavity.V(:,iCavity); % cavity volumes (state variable)
% arterial and venous length are equal
Len  = ArtVen.Len([1,1],:); Len=Len(:)'; % repesentative length of blood vessels
p0   = ArtVen.p0(:)'; % working pressure
A0   = ArtVen.A0(:)';

A      = max(1e-10,bsxfun(@rdivide,V,Len)); % vessel cross-section, safety A>0
ANorm  = bsxfun(@rdivide,A,A0); % cross-section normalized to wall volume

% fraction ap determines anti-collapse stiffness for small volume
ap    = 0.02; % with volume collapse to ap fraction -> steep slope pTrans<0
mk    = 1 + (P.ArtVen.k(:)'/3-2)/(1+ap); % stiffness exponential
pp    = (1+ap) * bsxfun(@power,ANorm,mk); % pressure with normal and large volume
pm    = -ap./ANorm; % anti-collapse pressure component at small volume
pTrans= bsxfun(@times, pp + pm, p0); % total pTrans
c0    = sqrt(0.25+ bsxfun(@times, bsxfun(@times,pp,mk)-pm, p0) / rhob );
%       c0=wave velocity with q=0, minimized to 0.5 m/s
Z0    = rhob * c0 ./ A; % wave impedance with flow=0

P.Cavity.pTrans(:,iCavity)= pTrans; % transmural pressure
P.Cavity.Z(:,iCavity)     = Z0    ; % wave impedance
P.Cavity.A(:,iCavity)     = A     ; % cross-section

end

