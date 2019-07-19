function pNodeVDot
% function pNodeVDot
% Pressures in cavities p -> pressure in Nodes: p
%                         -> flows to cavities: VDot
% Includes collapsibility of tubes
% Theo Arts, Maastricht University, Aug 7, 2013

global P

% flow of valves to/from nodes, initialization Node.q/Y/A
nTube  =P.Tube.n;

%=== index conversion sparse matrices
Valve2NodeProx= P.Valve.Valve2NodeProx;
Valve2NodeDist= P.Valve.Valve2NodeDist;
Cavity2Node   = P.Cavity.Cavity2Node;
Tube2NodeProx = P.Tube.Tube2NodeProx;
Tube2NodeDist = P.Tube.Tube2NodeDist;


%=== Valve flow contribution to node pressure
DqValve= P.Valve.q*(-Valve2NodeProx+Valve2NodeDist); %Added flow into node
DAValve= repmat(max(P.Valve.AOpen,P.Valve.ALeak) * ...
    (Valve2NodeProx+Valve2NodeDist),size(P.t)); %Added flow area around node

%==== AV flow through ArtVen
p0Av    = P.ArtVen.p0AV;    % reference AV pressure drop
q0Av    = P.ArtVen.q0AV;    % reference AV flow
kAv     = P.ArtVen.kAV;     % exponent non-linearity AV 'resistance'
iAv2CavP= P.ArtVen.iCavity; % ArtVen ->index-> proximal cavities
iAv2CavD= iAv2CavP+1;       % ArtVen ->index-> distal cavities
iAv2NdP = P.ArtVen.iNode;   % ArtVen -> index -> arterial Node
iAv2NdD = iAv2NdP+1;        % ArtVen -> index -> venous Node
pCav    = P.Cavity.p;       % Cavity pressures
pAvCavP = pCav(:,iAv2CavP); % ArtVen: proximal cavity pressure
pAvCavD = pCav(:,iAv2CavD); % ArtVen: distal   cavity pressure
Dp      = pAvCavP-pAvCavD;  % ArtVen: AV-pressure drop between Cavities

% Correction of ArtVen-flow resistance to volume change of 
% Artven-related cavities, Resistance ~ 1/Volume^2
DpV20     = p0Av.*P.ArtVen.Len.^2./sum(P.ArtVen.A0.^-2);
VCav      = P.Cavity.V;
DpV2      = abs(Dp)./(VCav(:,iAv2CavP).^-2 + VCav(:,iAv2CavD).^-2);
qNorm     = bsxfun(@power,bsxfun(@rdivide,DpV2,DpV20),kAv).*sign(Dp);
qAv       = bsxfun(@times,qNorm,q0Av); % AV flow
P.ArtVen.q= qAv;
Facq      = P.General.FacpControl; % FacpControl= current p / pRef

%===== Cavity contributions to node pressure and flow
iAv2Nd =[iAv2NdP,iAv2NdD];   % ArtVen -> index -> Node P+D
YCav  = 1./P.Cavity.Z;       % conductivity per cavity
qCav  = pCav .* YCav;        % pressure x conductivity per cavity
DYCav = YCav  * Cavity2Node; % internal Y conductance of flow source
DqCav = qCav  * Cavity2Node; % short circuit flow of node
DqCav(:,iAv2Nd)= DqCav(:,iAv2Nd) + [-qAv*Facq,qAv/Facq];
%      ArtVen flow + volume control
DACav = P.Cavity.A* Cavity2Node; % flow cross-sectional area

%===== Tube contributions to node pressure and flow
YTbProx= 1./P.Tube.ZR; % conductivity per tube
YTbDist= 1./P.Tube.ZL; % conductivity per tube
pTbProx= P.Tube.pProx; % source pressure
pTbDist= P.Tube.pDist; % source pressure
YTbS   =[YTbProx,YTbDist] ; % source conductivity
pTbS   =[pTbProx,pTbDist] ; % source pressure
qTbS   = pTbS .* YTbS;      % source flow to node
Tb2Nd  =[Tube2NodeProx;Tube2NodeDist]; % Tube -> matrix -> Node P+D
DYTube = YTbS * Tb2Nd; % added Tube related conductivity
DqTube = qTbS * Tb2Nd; % added Tube related source flow
DATube = P.Tube.A * ( Tube2NodeProx + Tube2NodeDist ); % added area

% No waterfall
YNode=           DYCav + DYTube; % total node conductivity
qNode= DqValve + DqCav + DqTube; % total node inflow with pNode=0
ANode= DAValve + DACav + DATube; % total area, seen from node

pNode= qNode ./ YNode; % node pressure by solving NodeInflow Dq=0

% Detect waterfall conditions in Tube
iTb2NdP= P.Tube.iNodeProx;
iTb2NdD= P.Tube.iNodeDist;
pTbExt = P.Tube.p - P.Tube.pTrans; % external tube pressure
iTb2Nd =[iTb2NdP,iTb2NdD];
pTbNd  = pNode(:,iTb2Nd);
YTbNd  = YNode(:,iTb2Nd);
pTbE   =[pTbExt,pTbExt];

% Waterfall by pressure drop
dpTb = dpWf(pTbNd,pTbE,pTbS,YTbS,YTbNd); %Waterfall pressure drop
DqTb = -YTbS.*dpTb;
pTbS = pTbS-dpTb; % waterfall
qNode= qNode + DqTb * Tb2Nd;

YTbProx= YTbS(:, 1:nTube       );
YTbDist= YTbS(:,(1:nTube)+nTube);
pTbProx= pTbS(:, 1:nTube       );
pTbDist= pTbS(:,(1:nTube)+nTube);

%Waterfall in Tube
pNode  = qNode ./ YNode; % Waterfall corrected node pressures

% d/dt Cavity volume
iCav2Nd      = P.Cavity.iNode; % iNode pointed by cavity
P.Cavity.VDot= (pNode(:,iCav2Nd)-P.Cavity.p).*YCav;

P.ArtVen.qProx=  P.Cavity.VDot(:,iAv2CavP)+qAv;
P.ArtVen.qDist= -P.Cavity.VDot(:,iAv2CavD)+qAv;

% Tube volume change
pNodeProx= pNode(:,iTb2NdP); %pressure Prox-node
pNodeDist= pNode(:,iTb2NdD); %pressure Dist-node
qP=(pNodeProx-pTbProx).*YTbProx;
qD=(pTbDist-pNodeDist).*YTbDist;
P.Tube.qProx= qP;
P.Tube.qDist= qD;
P.Tube.VDot = qP-qD;

P.Node.p = pNode;
P.Node.Y = YNode;
P.Node.q = qNode;
P.Node.A = ANode;

end

function dp=dpWf(pN,pE,pS,YS,YN)
% {pN,pE,pS}= pressure {Node, External, Source}
% {Y,YN}= conductivity {waterfall, Node total}
% dp= decrease of source pressure to simulate waterfall
%=====================
eps= 0.2;
PS = pS-pE;
PN = pN-pE;
OutFlow= double(PS>PN); % With inflow no waterfall
% p  = OutFlow.*max(0,min(0,PS)-PN); %Waterfall
p  = OutFlow.*max(0,min(0,PS)-PN); %Waterfall
x  = 1-YS./YN; % effect of node impedance <-> tube impedance
z = x+eps*exp(-(x/eps).^2); %safety to avoid zero division
dp= p./z;
end

