function TubeV2p
% function TubeV2p
% Tube volume V -> transmural pressure pTrans, wave propagation velocity
% c0, proximal and distal source pressure pProx, pDist and source
% impedances ZR (prox), ZL (dist), using delayed pressures.
% Theo Arts, Maastricht University, Oct 13, 2012

% Imposed changes in the calculation of vicous friction of blood flow, now based on an
% approximation that is in agreement with Womersley theory, Maarten Heusinkveld, August 2018

global P

Len   = P.Tube.Len   ; % repesentative length of blood vessels
p0    = P.Tube.p0    ; % working pressure
A0    = P.Tube.A0    ;
rhob  = P.General.rhob;

A      = max(1e-10,bsxfun(@rdivide,P.Tube.V,Len)); % vessel cross-section
ANorm  = bsxfun(@rdivide,A,A0); % cross-section normalized to physiologic volume

eta    = 0.0030;
hRPois = bsxfun(@rdivide,4*pi*eta*Len, A.^2);

% calculation [pTrans, c0, Z0] with anti-collaps stiffness
ap    = 0.02; % introduces: Small volume -> negative transmural pressure pTrans
mk    = 1 + (P.Tube.k/3-2)/(1+ap); % stiffness exponential
pp    = (1+ap) * bsxfun(@power, ANorm, mk);
pm    = -ap./ANorm; % anti-collapse pressure
pTrans= bsxfun(@times,pp + pm, p0);
c0    = sqrt( bsxfun(@times, bsxfun(@times,pp,mk)-pm, p0) / rhob );
Z0    = rhob * c0 ./ A; % wave impedance with tube flow=0

if strcmp(P.General.Profile,'Approximate')
    [Att,cWom,Z0] = CalcApproximateVelocityProfile(c0,'Tube');   
elseif strcmp(P.General.Profile,'Flat')
    Att           = zeros(length(P.t),length(Len));%
    cWom          = c0;
end

% Flow dependency of wave velocity and impedance
% getting delayed signals for source pressure pProx and pDist
Dt= P.General.Dt;
it= round((P.t-P.General.tStart)/Dt)+1;

q     = P.Tube.q(it,:); % tube flow
ZR = Z0; ZL = Z0; cR = cWom; cL = cWom;
% Calculate the pressure drop due to Poiseuille friction 
DpPois = q.*hRPois; % effect DC

% Apply correction for Poiseuille resistance effect      
P.Tube.pProx  = P.Tube.pP(it,:)+DpPois;
P.Tube.pDist  = P.Tube.pD(it,:)-DpPois;
P.Tube.ZR     = ZR     ; %R-wave impedance % Note in this case assumed equal
P.Tube.ZL     = ZL     ; %L-wave impedance % Note in this case assumed equal
P.Tube.cL     = cL     ; %R-wave velocity  % Note in this case assumed equal
P.Tube.cR     = cR     ; %L-wave velocity  % Note in this case assumed equal
P.Tube.A      = A      ;
P.Tube.pTrans = pTrans ;
P.Tube.p      = pTrans ;
P.Tube.Att    = Att    ; % wave attenuation constant per unit length [1/m]

end

