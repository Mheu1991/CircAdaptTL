function TubeV2p
% function TubeV2p
% Tube volume V -> transmural pressure pTrans, wave propagation velocity
% c0, proximal and distal source pressure pProx, pDist and source
% impedances ZR (prox), ZL (dist), using delayed pressures.
% Theo Arts, Maastricht University, Oct 13, 2012

% Imposed changes in the calculation of vicous friction of blood flow, now based on an
% approximation that is in agreement with Womersley theory (considered gold
% standard), Maarten Heusinkveld, August 2018

global P

Len   = P.Tube.Len   ; % repesentative length of blood vessels
p0    = P.Tube.p0    ; % working pressure
A0    = P.Tube.A0    ;
rhob  = P.General.rhob;

A      = max(1e-10,bsxfun(@rdivide,P.Tube.V,Len)); % vessel cross-section
ANorm  = bsxfun(@rdivide,A,A0); % cross-section normalized to physiologic volume

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
vb    = q./A; % mean blood velocity
hvdc0 = 0.5*vb./cWom; % blood velocity, normalized to wave velocity
% Simple general derivation
b = sqrt(1+hvdc0.^2);
bR= b+hvdc0; bL= b-hvdc0;
cR= sqrt(0.25+(cWom.*bR).^2);% wave velocity clipped to >0.5 m/s
cL= sqrt(0.25+(cWom.*bL).^2);% wave velocity clipped to >0.5 m/s
ZR= Z0./bR ; ZL= Z0./bL;

P.Tube.pProx  = P.Tube.pP(it,:);
P.Tube.pDist  = P.Tube.pD(it,:);
P.Tube.ZR     = ZR     ; %R-wave impedance
P.Tube.ZL     = ZL     ; %L-wave impedance
P.Tube.cL     = cL     ; %R-wave velocity
P.Tube.cR     = cR     ; %L-wave velocity
P.Tube.A      = A      ;
P.Tube.pTrans = pTrans ;
P.Tube.p      = pTrans ;
P.Tube.Att    = Att    ; % wave attenuation constant per unit length [1/m]

end

