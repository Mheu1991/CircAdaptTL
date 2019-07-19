function TimingInit
% function Timing
% Contains definition of depolarization pathways
% Sets patch depolarization times for 1 beat
% Theo Arts, Maastricht University, April 30, 2013

global P
G        = P.General    ;
tCycleRef= G.tCycleRest ; % Reference tCycle for body size scaling
tCycle   = G.tCycle     ; % current cycle time
TimeFac  = G.TimeFac    ; % scaling of contraction time
TauAv    = G.TauAv      ; % AV-delay, controlled in Adapt0

%=== setting TimeAct [s] duration of contraction for Ls==LsStress0Act
% in atrial and ventricular wall patches
ta= 0.15*(tCycle/0.85)*TimeFac; % atrial activation duration
tv= (0.10*tCycleRef +0.40*tCycle)*TimeFac; % ventricular " "
% Find atrial patches, discrimination atria<->ventricles
iPa=[]; iPv=[];
for i=1:P.Patch.n
    if regexp(P.Patch.Name{i},'a')==2
        iPa=[iPa,i];
    else
        iPv=[iPv,i];
    end
end
P.Patch.TimeAct(iPa)=ta;
P.Patch.TimeAct(iPv)=tv;
%================ end setting TimeAct====================

%====================================================
% definition of interpatch conduction pathways
P.Net.Depolarization.Pathway={...
    'Ra1','Ra1',tCycle;...
    'Ra1','La1',0.02 * (tCycleRef/0.85) * TimeFac;...
    'Ra1','Rv1',TauAv;...
    'Ra1','Lv1',TauAv;...
    'Rv1','Sv1',0;...
    'Sv1','Lv1',0;...
    };
P.Net.Depolarization.Pace         = 'Ra1'; % leading pacemaker patch
P.Net.Depolarization.TauInterPatch= 0.005; % default interpatch delay
P.Net.Depolarization.TauRefrac    = 0.25 ; % default refractory period

% === CONDUCTION PATHWAYS as defined in P.Net.Depolarization
Pathway= P.Net.Depolarization.Pathway;
nP     = length(Pathway); %number of paths
DepPath= zeros(nP,3); % initialization of [iPatchProx,iPatchDist,Tau]
for i=1:nP; % depolarization pathways copied from P.Net.Depolarization
    Path        = Pathway(i,:);
    DepPath(i,:)= [Get('Patch','Index',Path(1:2)),Path{3}];
end
% Intrawall interpatch pathways
DepP    = [];
Tau     = P.Net.Depolarization.TauInterPatch; %Default interpatch delay
for iW  = 1:P.Wall.n
    iP  = P.Wall.iPatch(iW);
    nP  = P.Wall.nPatch(iW);
    Aux = (0:nP-2)';
    dM  = [iP+Aux,(iP+1)+Aux];
    DepP=[DepP;dM;fliplr(dM)];
end
P.Patch.DepPath=[ DepPath; [DepP,repmat(Tau,[size(DepP,1),1])] ]';
% indexed depolarization pathways
P.Patch.TauRefrac=repmat(P.Net.Depolarization.TauRefrac,[1,P.Patch.n]);
% default refractory period
P.Patch.iPace= Get('Patch','Index',P.Net.Depolarization.Pace);
% Pacemaker patch
% === END indexation CONDUCTION PATHWAYS ================

end

