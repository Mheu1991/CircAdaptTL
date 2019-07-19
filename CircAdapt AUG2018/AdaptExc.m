function AdaptExc
%function AdaptExc;
% Specific adaptation to exercise conditions
% Control of:
% Blood vessel wall thickness by wall stress
% Heart wall thickness and cavity volume by myocardial tissue load
% Bag pressure
% Theo Arts, Maastricht University, Feb 20, 2013

global P

P.General.In=[P.General.In; P2StatVec('In')];
%=== Begin adapt special
%=== Actions of adaptation at excercise
%=== Adapt ArtVen wall thickness and Patches
ArtVenAdapt('All','WallVolume');
TubeAdapt('All','WallVolume');
PatchAdapt('All','All');

% % %==== Fixation Septal/Lv wall area to 1/2
% since adaptation is critical at that point +++++++
iP=Get('Wall','iPatch',{'Lv','Sv'});
if ~isempty(iP)
    nP=Get('Wall','nPatch',{'Lv','Sv'});
    AmRefLv=P.Patch.AmRef(iP(1)-1+(1:nP(1)));
    AmRefSv=P.Patch.AmRef(iP(2)-1+(1:nP(2)));
    LvDSv=2;
    Fac=sqrt((sum(AmRefLv)/sum(AmRefSv))/LvDSv);
    P.Patch.AmRef(iP(1)-1+(1:nP(1)))=AmRefLv/Fac;
    P.Patch.AmRef(iP(2)-1+(1:nP(2)))=AmRefSv*Fac;
    disp(['VWall Lv+Sv= ',...
        num2str(sum(Get('Wall','VWall',{'Lv','Sv'}))*1e6,'%10.1f')]);
end
% % %---- End fixation (AmRefLv/AmRefSv)->2 -----

BagAdapt;
%=== end actions of adaptation


P.General.Out=[P.General.Out;P2StatVec('Out')];
if P.General.Fast;% Put new start values in structure P
    StatVec2P(SteadyStateP);
end


end

%==== Specific adaptation functions ==========================

%===== Stationarity functions for Exercise =====+++++++++++++++
function Vec= P2StatVec(InOut) % vector of variables used for reaching stationarity
global P;
VecV=Get('Cavity','V','All');
if strcmp(InOut,'In')
    Vec=VecV(1,:);
else
    Vec=VecV(end,:);
end
end

function StatVec2P(Vec) % Inverse of P2StatPar
% fills values at the end
global P;
nC=P.Cavity.n;
Put('Cavity','V','All',Vec(1:nC));
% P.Patch.SfPas=Vec(nC+1:end);
end

function BagAdapt %  Bag adaptation to pAdapt reference
global P
if P.Bag.n==0; return; end
pMax = max(P.Bag.pTrans);
Fac_pMax= max(0.1,pMax./P.Bag.pAdapt);
P.Bag.VRef= P.Bag.VRef.*(Fac_pMax.^(0.3 ./ P.Bag.k));
disp(P.Bag.Name);disp(log(Fac_pMax));
end

