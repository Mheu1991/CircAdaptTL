function pCavity
% function pCavity
% Bags enclose Cavities, Walls, Tubes and other Bags
% Transmural pressures are added to absolute pressures
% in Cavities, Tubes and other Bags
% Bags may be nested
% Theo Arts, Maastricht University, Oct 13, 2012

global P
nt=length(P.t);

if P.Bag.n>0
    VC= P.Bag.VC;
    VW= P.Bag.VW;
    VT= P.Bag.VT;
    VB= P.Bag.VB;
    P.Bag.V= P.Cavity.V*VC+P.Tube.V*VT+repmat(P.Wall.VWall*VW,[nt,1]);
    
    % Bag Volume -> Bag pTrans
    VRef        = P.Bag.VRef  ;
    pAdapt      = P.Bag.pAdapt;
    P.Bag.pTrans= bsxfun(@times,...
        bsxfun(@power, bsxfun(@rdivide,P.Bag.V,VRef),P.Bag.k),...
        pAdapt);

    % Pressure increments due to Bag-pTrans
    % These increments are external pressures on Cavity, Tube and Bag
    dPC= P.Bag.pTrans * VC'; %Cavities
    dPT= P.Bag.pTrans * VT'; %Tubes
    dPB= P.Bag.pTrans * VB'; %Bags

    % Bags
    P.Bag.p= P.Bag.pTrans + dPB;
else
    dPC = 0; dPT = 0;
end

%Cavities
%Default no external pressure
P.Cavity.p= P.Cavity.pTrans + dPC; %Default no external pressure

%Tubes
%Default no external pressure
if P.Tube.n>0 % Add dPT to mean tube and source pressures
    P.Tube.p    =P.Tube.pTrans + dPT;
    P.Tube.pProx=P.Tube.pProx  + dPT;
    P.Tube.pDist=P.Tube.pDist  + dPT;
end

% initialization node pressure calculation
nNode   = P.Node.n;
Aux     = zeros(nt,nNode);
P.Node.p= Aux;
P.Node.Y= Aux;
P.Node.q= Aux;
P.Node.A= Aux;

end

