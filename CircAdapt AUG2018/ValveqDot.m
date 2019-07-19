function ValveqDot
% function ValveqDot
% Node pressure differences -> valve flow acceleration qDot
% Valves AvValves refer to ventricular valves with papillary muscles
% Modification to more robustness of valve closure
% Theo Arts, Maastricht University, Mar 21, 2016

global P

Dt       = P.General.Dt;
iNodeProx= P.Valve.iNodeProx;
iNodeDist= P.Valve.iNodeDist;
q        = P.Valve.q  ; % flow
nt       = size(q,1);
AOpen    = repmat(P.Valve.AOpen,[nt,1]); % open valve cross-section
ALeak    = repmat(P.Valve.ALeak,[nt,1]); % closed valve cross-section
Len      = P.Valve.Len; % effective length of flow channel
rhob     = 1050       ; % density of blood
%=================
% allows AV-valve diastolic regurgitation
Ws   = P.Valve.AvWalls; Vs=P.Valve.AvValves;
T    = P.Wall.T(:,Ws); %wall tension related to pap. muscle
DADT = P.Wall.DADT(:,Ws); %wall stiffness
Am0  = P.Wall.Am0(:,Ws); %wall zero-stress area
Diast= tanh( 30*max(0,T.*DADT./Am0-0.1).^2 ); % diastole->1.0
ALeak(:,Vs)= 0.3*(AOpen(:,Vs)-ALeak(:,Vs)).*Diast+ALeak(:,Vs); 
%    Large AV-valve leak in diastole

Dp    = P.Node.p(:,iNodeProx)-P.Node.p(:,iNodeDist); % pressure drop

AMax  = max(AOpen,ALeak);
AProx = P.Node.A(:,iNodeProx)-AMax;
ADist = P.Node.A(:,iNodeDist)-AMax;

% Slow closure to avoid numerical problems
Bwq = +(q<0); % Backward flow
Bwp = +(Dp<0); % Backward pressure
Fw  = 1-Bwq.*Bwp; % Forward valve leaflet opening
A   = ALeak + Fw.*(AOpen-ALeak); % Flow cross-section

%   Bernouilli pressure drop over valve orifice
R  = 1.0+rhob*abs(q).*abs( A.^-2 - 1./((1-Bwq).*AProx.^2+Bwq.*ADist.^2) );
DpB= 0.5 .* q .* R; % Bernouilli pressure drop
L  = 1.5*rhob*( bsxfun(@ldivide,A,Len) ...
    + 0.5*(1./sqrt(AProx)+1./sqrt(ADist))); % inertia
Tau = L./R; % time constant of qDot, used to avoid numerical problems 
% related to short time constant of valve closure
% qDot= (1-exp(-Dt./Tau)).*(Dp-DpB)./(Dt*R); %linear resistance
qDot= (1-1./(1+Dt./Tau)).*(Dp-DpB)./(Dt*R); %non-lin Bernouilli solution

P.Valve.qDot= qDot; % flow derivative
end

