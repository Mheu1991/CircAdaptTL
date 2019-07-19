function OutDotT=SVarDot(tDummy,SVarT,flag)
% function OutDotT=SVarDot(tDummy,SVarT,flag)
% tDummy= time (not in use)
% SVarT = column vector of state variables SVar (=row vector)
% flag  = not in use
% OutDotT= column vector of derivatives
% State Variables -> Time derivatives of State Variables
% Ready to be used in MatLab function odeXXX() to solve set of Diff Eq
% Theo Arts, Maastricht University, Oct 13, 2012
%====

global P

P.SVar= SVarT'; % store state variables SVar
SVar2P; % state variables SVar -> physiologic representation P.xx
P.tDot=ones(size(P.t)); % time derivative of time = 1
MemAlloc; % necessary memory allocations
ArtVenV2p; % arteries to veins compliant network
PatchWallA2T; % patch and wall: Am= Am0+T*DADT
ChamberV2p; % Chamber, pTrans and Wall.Am,T
TriSegV2p; % TriSeg, pTrans and Wall.Am,T
Wall2Patch2Sarc; % filling Sarc with LsiDot,CDot
TubeV2p; % Delayed source transmural pressures and resistances of tube
pCavity; % Determines p in cavities by adding bag pressures
pNodeVDot; % Node pressures and cavity VDot's
ValveqDot; % flow derivatives qDot
TubeDelays; % needed for Tube-delays

P2SVarDot; % transfer of derivatives to P-structure

OutDotT= real(P.SVarDot)'; % odeXX requires OutDotT to be a column vector

end

