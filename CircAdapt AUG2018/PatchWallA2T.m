function PatchWallA2T
% function PatchWallA2T
% Determines linear function patch area Am = fu Tension T
% Am(T)=Am0+T*DADT
% Constants Am0 and DADT are determined from Patch state variables
% and reference geometry
% Theo Arts, Maastricht University, Aug 31, 2014

global P;

% Patch Am= Am0+DADT*T, provides Am0 and DADT
AmRef= P.Patch.AmRef; % midwall area for ref sarcomere length 2mu
LsRef= P.Patch.LsRef; % ref sarcomere length 2mu
VWall= P.Patch.VWall; % wall volume

Lsi   = P.Patch.Lsi; % unloaded sarc length= state variable
Lambda= bsxfun(@rdivide,Lsi      ,LsRef); %extension factor
Am    = bsxfun(@times  ,Lambda.^2,AmRef); %actual midwall area
Ef    = log(Lambda); %natural fiber strain

P.Patch.Ef  = Ef; % fiber strain Ef with zero length SE-element

SarcEf2Sf; % sarcomere strain->stress

Sf    = P.Patch.Sf; % sarcomere stress
DSfDEf= P.Patch.DSfDEf;

%Aux   = DSfDEf - 2*Sf; % Theoretical value, not robust for low stiffness
Aux1   = 2*Sf./DSfDEf;
Aux    = DSfDEf.*(1-Aux1./sqrt(1+Aux1.^2));
% Approximates theoretical value of Aux, but zero division is avoided

DADT  = 4 * Am.^2 ./ bsxfun(@times,Aux,VWall);
Am0   = Am - 2*Am.*Sf./Aux;
P.Patch.DADT= DADT; % area compliance
P.Patch.Am0 = Am0; % zero tension area
% P.Patch.T   = bsxfun(@times, Sf./Am, (0.5 *VWall)); % tension

% Wall is composed of patches: Also for wall: Am(T)=Am0+DADT*T
for iWall=1:P.Wall.n
    iPatch= (P.Wall.iPatch(iWall)-1)+(1:P.Wall.nPatch(iWall));
    P.Wall.VWall(iWall)  = sum(P.Patch.VWall(iPatch));
    P.Wall.Am0(:,iWall)  = P.Wall.AmDead(iWall)+sum(P.Patch.Am0(:,iPatch),2);
    P.Wall.DADT(:,iWall) = sum(P.Patch.DADT(:,iPatch),2);
end
end

