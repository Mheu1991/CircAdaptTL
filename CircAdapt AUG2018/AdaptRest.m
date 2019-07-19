function AdaptRest
% function AdaptRest
% Adaptation at rest conditions
% Adaptation of:
%   Vessel cross-section
%   If needed valve diameter
% Theo Arts, Maastricht University, April 2, 2014

global P
ValveDiamAdapt=0; %==0/1: Valve orifice adaptation No/Yes

P.General.In=[P.General.In; P2StatVec('In')];
%=== Begin adapt special
ArtVenAdapt('All','Diameter');
TubeAdapt(  'All','Diameter');

if ValveDiamAdapt;% Adapt diameter of valves to the connected blood vessel
    Put('Valve','AOpen',{'VcRa','RvPuAr','PuVeLa','LvAo'},...
        mean(Get('Node','A',{'Vc','PuAr','PuVe','Ao'}))-...
        Get('Valve','AOpen',{'VcRa','RvPuAr','PuVeLa','LvAo'}));
    Put('Valve','ALeak',{'VcRa','PuVeLa'},...
        mean(Get('Valve','AOpen',{'VcRa','PuVeLa'})) );
    % mitral and tricuspid valve are larger
    Put('Valve','AOpen',{'RaRv','LaLv'},...
        1.5*Get('Valve','AOpen',{'RvPuAr','LvAo'}) );
end
P.General.Out=[P.General.Out;P2StatVec('Out')];
if P.General.Fast;% Put new start values in structure P
    StatVec2P(SteadyStateP);
end
end


%==== Specific adaptation functions ==========================


%===== Stationarity functions for AdaptRest =====
function Vec= P2StatVec(InOut) % vector of variables used for reaching stationarity
VecV=Get('Cavity','V','All');
if strcmp(InOut,'In')
    Vec= VecV(1,:);
else
    Vec= VecV(end,:);
end
end

function StatVec2P(Vec) % Inverse of P2StatPar
% fills values at the end
global P;
nC=P.Cavity.n;
Put('Cavity','V','All',Vec(1:nC));
end

