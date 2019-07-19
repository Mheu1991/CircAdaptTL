function LinkMatrix
% function LinkMatrix
% Links between elements are stored in matrices
% Matrices are constructed, based on information, described by strings
% Main purpose: matrices speeds up calculation, because string reading is slow
% Theo Arts, Maastricht University, Oct 13, 2012

global P

%sparse matrix for Cavity <-> Node connections
P.Cavity.Cavity2Node=sparse(zeros(P.Cavity.n,P.Node.n));
Aux=(1:P.Cavity.n)+(P.Cavity.iNode-1)*P.Cavity.n; % index
P.Cavity.Cavity2Node(Aux)=1;

% sparse matrix for ArtVen <-> Node connections
Aux= sparse(zeros(P.ArtVen.n,P.Node.n));
P.ArtVen.ArtVen2NodeArt=Aux;
P.ArtVen.ArtVen2NodeVen=Aux;
AuxArt=(1:P.ArtVen.n)+(P.Cavity.iNode(P.ArtVen.iCavity  )-1)*P.ArtVen.n;
AuxVen=(1:P.ArtVen.n)+(P.Cavity.iNode(P.ArtVen.iCavity+1)-1)*P.ArtVen.n;
P.ArtVen.ArtVen2NodeArt(AuxArt)=+1;
P.ArtVen.ArtVen2NodeVen(AuxVen)=+1;

% sparse matrix for Valve <-> Node connections
Aux= sparse(zeros(P.Valve.n,P.Node.n));
P.Valve.Valve2NodeProx=Aux;
P.Valve.Valve2NodeDist=Aux;
AuxProx=(1:P.Valve.n)+(P.Valve.iNodeProx-1)*P.Valve.n;
AuxDist=(1:P.Valve.n)+(P.Valve.iNodeDist-1)*P.Valve.n;
P.Valve.Valve2NodeProx(AuxProx)=+1;
P.Valve.Valve2NodeDist(AuxDist)=+1;

% sparse matrix for Tube <-> Node connections
Aux= sparse(zeros(P.Tube.n,P.Node.n));
P.Tube.Tube2NodeProx=Aux;
P.Tube.Tube2NodeDist=Aux;
AuxProx=(1:P.Tube.n)+(P.Tube.iNodeProx-1)*P.Tube.n;
AuxDist=(1:P.Tube.n)+(P.Tube.iNodeDist-1)*P.Tube.n;
P.Tube.Tube2NodeProx(AuxProx)=+1;
P.Tube.Tube2NodeDist(AuxDist)=+1;

% Valve-Wall connections (papillary muscles)
P.Valve.AvValves=Get('Valve','Index',{'LaLv','RaRv'});
P.Valve.AvWalls =Get('Wall','Index',{'Lv','Rv'});

% Bag indexings
% sparse matrices for volume contributions to Bag by:
% Cavity-VC, Wall-VW, Tube-VT concatenated to VCWT
% sparse matrices for pressure contributions from Bag to:
% Cavity-PC, Wall-PW, Tube-PT
nBag= P.Bag.n;
VC=zeros(P.Cavity.n,nBag);
VW=zeros(P.Wall.n  ,nBag);
VT=zeros(P.Tube.n  ,nBag);
VB=zeros(P.Bag.n   ,nBag);
for iBag=1:nBag
    VC(P.Bag.iCavity{iBag},iBag)=1; %Cavity volume index
    VW(P.Bag.iWall{iBag}  ,iBag)=1; %Wall volume index
    VT(P.Bag.iTube{iBag}  ,iBag)=1; %Tube volume index
    VB(P.Bag.iBag{iBag}   ,iBag)=1; %Bag volume index
end
P.Bag.VC= VC;
P.Bag.VW= VW;
P.Bag.VT= VT;
P.Bag.VB= VB;

%==== Prepare handling of waves in tubes
Dt    = P.General.Dt;
nt    = round(P.General.tCycle/Dt);
P.General.tCycle= nt*Dt; % set integer number of time steps per cycle
nTube = P.Tube.n;
Mat   = zeros(nt,nTube);
ntOld = size(P.Tube.pL,1);
ntm1  = min(nt,ntOld)-1; % checks on new tCycle, modifies delay matrices
pL= Mat; % left wave
pR= Mat; % right wave
pL(end-ntm1:end,:)= P.Tube.pL(end-ntm1:end,:);
pR(end-ntm1:end,:)= P.Tube.pR(end-ntm1:end,:);
P.Tube.pL= pL; % circular delay memory
P.Tube.pR= pR; % circular delay memory
P.Tube.pP= P.Tube.pP([end,end],:); % proximal source pressures
P.Tube.pD= P.Tube.pD([end,end],:); % distal source pressures
P.Tube.q = P.Tube.q([end,end],:);
end

