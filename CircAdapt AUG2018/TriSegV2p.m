function TriSegV2p
% function TriSegV2p
% TriSeg is a 3-wall structure (Left,Septal,Right) with 2 cavities (R,L)
% Calculates: cavity volumes V -> dimensions of the 'double bubble',
% myofiber stress Sf, wall tension T and cavity pressures p
% VS and YS repesent septal volume displacement and junction radius.
% State variables P.TriSeg V and Y represent initial estimates to calculate
% VS and YS accurately.
% Theo Arts, Maastricht University, Oct 13, 2012

global P;
if P.TriSeg.n==0 %if there is no chamber
    return
end

TriSeg=P.TriSeg;

n  = TriSeg.n; % number of TriSeg's
iCavity = TriSeg.iCavity; %related cavities
iWall   = TriSeg.iWall  ; %related walls
rhob    = P.General.rhob; % blood density

for i=1:n %for all TriSeg's
    iC   = iCavity(:,i)+(0:1); % 2 cavities
    iW   = iWall(:,i)  +(0:2); % 3 walls
    Am0  = P.Wall.Am0(:,iW); %zero stress wall area
    DADT = P.Wall.DADT(:,iW); %wall compliance
    nt   = size(Am0,1); %number of time points
    VWall= P.Wall.VWall(iW); % 3 wall volumes
    VWL  = mean(VWall(1:2)); % enclosed wall volume 1st cavity
    VWR  = mean(VWall(2:3)); % enclosed wall volume 2nd cavity
    VCav = P.Cavity.V(:,iC);
    V    = max(0,bsxfun(@plus,VCav,[VWL,VWR])); %2 midwall volumes(t)

    VS= TriSeg.V(:,i);
    YS= TriSeg.Y(:,i); % 1st estimate of [V,Y]
    YRef= mean(YS); VRef= YRef^3;
    dvR = 0.02;     dyR = dvR;
    dv=dvR*VRef;
    dy=dyR*YRef; % increments for d[Txy]/dVY

    as=sin(pi/12); ac=cos(pi/12); bs=sin(pi/4);
    % coefficients to for a equilateral triangle in 2D
    [Tx0a,Ty0a]=VY2Txy(VS-bs*dv,YS+bs*dy,Am0,DADT,V); % tension Txy =fu(VY)
    [TxVa,TyVa]=VY2Txy(VS+ac*dv,YS+as*dy,Am0,DADT,V); % partial V derivative
    [TxYa,TyYa]=VY2Txy(VS-as*dv,YS-ac*dy,Am0,DADT,V); % partial Y derivative
    % matrix of coefficxients to determine partial derivatives
    DX=[-bs  bs
         ac  as
        -as -ac];
    Ddvy= pinv(DX')*diag([1/dv,1/dy]);
    DTx= [Tx0a,TxVa,TxYa]*Ddvy; % [dTx/dx,dTx/dy]
    DTy= [Ty0a,TyVa,TyYa]*Ddvy; % [dTy/dx,dTy/dy]

    DET= DTx(:,1).*DTy(:,2)-DTx(:,2).*DTy(:,1); % determinant
    Tx0=mean([Tx0a,TxVa,TxYa],2);
    Ty0=mean([Ty0a,TyVa,TyYa],2);
    dV=(-DTy(:,2).*Tx0+DTx(:,2).*Ty0)./DET;
    dY=(+DTy(:,1).*Tx0-DTx(:,1).*Ty0)./DET;
    VS=VS+dV;
    YS=YS+dY;
    Tau=TriSeg.Tau;
    TriSeg.VDot(:,i)= dV/Tau; % keep track of solution
    TriSeg.YDot(:,i)= dY/Tau;

    for j=1:1 % extra iterations for solution of TriSeg geometry
        [Tx0,Ty0]=VY2Txy(VS,YS,Am0,DADT,V);
        dV=(-DTy(:,2).*Tx0+DTx(:,2).*Ty0)./DET;
        dY=(+DTy(:,1).*Tx0-DTx(:,1).*Ty0)./DET;
        VS=VS+dV;
        YS=YS+dY;
    end
    
    % writing geometric data TriSeg
    TriSeg.YS(:,i)=YS; % final solution Rv-Sv-Lv junction radius
    TriSeg.VS(:,i)=VS; % volume of septal cap
    
    % writing wall data
    [Tx0,Ty0,Am,Cm,T]=VY2Txy(VS,YS,Am0,DADT,V); % solution TriSeg
    P.Wall.Am(:,iW)= max(Am0,Am); % wall area
    P.Wall.Cm(:,iW)= Cm; % wall curvature
    P.Wall.T(:,iW) = max(0,T) ; % wall tension
    pTrans= 2*Cm.*T; % transmural pressure
    P.Wall.pTrans(:,iW)=pTrans;
     
    % Cavity impedance properties, needed to make node connection +++++++
    % efficienter maken
    Vw= repmat([VWL,VWR],[nt,1]);
    Vm= V + Vw;
    Len= 2*Vm.^(1/3);
    A  = ( V + 0.1*Vw ) ./Len;
    Z0 = sqrt(rhob ./ abs(A.*DADT(:,[1,3]).*Len)); %Cavity wave impedance
    P.Cavity.A(:,iC) = A; % cross-sectional area for valve inflow and outflow pressure
    P.Cavity.Z(:,iC) = Z0;
    P.Cavity.pTrans(:,iC)= [-pTrans(:,1),pTrans(:,3)];
    
end

P.TriSeg=TriSeg;
end
function [Tx,Ty,Am,Cm,Tm]=VY2Txy1(VS,YS,Am0,DADT,VLR)
% 1st order approximation of TriSeg according to J. Lumens et al.
% For a wall with zero-stress area Am0 and compliance DADT
% cap volume VS, junction radius YS and cavity volumes VLR=[VL,VR]
% Result: Summed axial and radial tension components [Tx,Ty] on junction
% for 3 walls: Am= midwall area, Cm= midwall curvature, Tm= wall tension
ARef=YS.^2;
VRef=YS.^3;
VS  = VS./VRef;
Am0 = bsxfun(@rdivide,Am0 ,ARef);
DADT= bsxfun(@rdivide,DADT,ARef);
VLR = bsxfun(@rdivide,VLR ,VRef);
Vm  = [VLR,VS] * [...
    -1     0     0
     0     0     1
     1     1     1];
% Solving 3rd order polynomial
% Mathematica: Solve[x^3 + 3x  - 2V == 0, x]
% Q= (V + Sqrt(V^2 + 1))^(1/3);  x= Q - 1/Q;
SignVm= sign(Vm); Vm=abs(Vm);
V     = (3/pi)*Vm;
Q     = (V + sqrt(V.^2 + 1)).^(1/3);
Xm    = SignVm .* ( Q - 1./Q );

%calculate midwall area Am and curvature Cm=1/rm
X2    = Xm.^2;
R2    = X2+1;
Am    = pi*R2; % midwall cap area, buckling with T<0%+++++
Cm    = 2*Xm./R2; % midwall cap curvature

% calculation of tension T and components Tx, Ty
Tm  = (Am-Am0)./DADT;
Sin = Cm;
Cos = (1-X2)./R2;
Txi = Cos.*Tm; %
Tyi = Sin.*Tm;
Tx  = sum(Txi,2); % normalized axial tension component
Ty  = sum(Tyi,2); % normalized radial tension component
Am  = bsxfun(@times  ,Am,ARef);
Cm  = bsxfun(@rdivide,Cm,YS  );
end

function [Tx,Ty,Am,Cm,Tm]=VY2Txy(VS,YS,Am0,DADT,VLR)
% 1st order approximation of TriSeg according to J. Lumens et al.
% For a wall with zero-stress area Am0 and compliance DADT
% cap volume VS, junction radius YS and cavity volumes VLR=[VL,VR]
% Result: Summed axial and radial tension components [Tx,Ty] on junction
% for 3 walls: Am= midwall area, Cm= midwall curvature, Tm= wall tension
Vm= [VLR,VS]* [...
    -1     0     0
     0     0     1
     1     1     1];
Ym=[YS,YS,YS];

% Solving 3rd order polynomial
% Mathematica: Solve[x^3 + 3y^2x  - 2V == 0, x]
% Q= (V + Sqrt(V^2 + y^6))^(1/3);  x= Q - y^2/Q;
SignVm= sign(Vm); Vm=abs(Vm);
V     = (3/pi)*Vm;
Q     = (V + sqrt(V.^2 + Ym.^6)).^(1/3);
Xm    = SignVm .* ( Q - Ym.^2 ./ Q );

%calculate midwall area Am and curvature Cm=1/rm
X2    = Xm.^2; Y2= Ym.^2;
R2    = X2+Y2;
% Am    = pi*R2; % midwall cap area
Am    = pi*R2; % midwall cap area, buckling with T<0
Cm    = 2*Xm./R2; % midwall cap curvature

% calculation of tension T and components Tx, Ty
Tm=(Am-Am0)./DADT;
Sin= Ym.*Cm;
Cos= (Y2-X2)./R2;
Txi = Cos.*Tm; %
Tyi = Sin.*Tm;
TRef=sqrt(sum(Tm.^2,2)); 
Tx= sum(Txi,2)./TRef; % axial tension component
Ty= sum(Tyi,2)./TRef; % radial tension component
end

