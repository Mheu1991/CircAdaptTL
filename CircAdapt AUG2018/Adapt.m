function Adapt
%function Adapt
% Common to all adaptation procedures.
% Embeds more dedicated adaptation procedures
%    executed between beats:
% Blood pressure is controlled by change of circulating volume
% Flow is set by P.General.SaturationControl:
%     0-> target flow per ArtVen-element
% or: 1-> flow is determined by oxygen usage= AV-Dsaturation x flow
% Theo Arts, Maastricht University, April 26, 2014

global P
save PTemp P; %saves last intermediate solution

% Systemic blood flow adjusted per ArtVen-element by resistance change
% Control of systemic pressure by adjustment of circulatory blood volume
FbFac=P.General.AdaptFeedback; % Feedback constant. If FbFac==0, no control

% Finding systemic and pulmonary ArtVen's
aux= regexpi(lower(P.ArtVen.Name),'pu'); %all names with 'pu' are pulmonary
Sy = P.ArtVen.Name(cellfun('isempty',aux));% {'Ca','Br','Ce','Fe'} systemic ArtVens
Pu = P.ArtVen.Name(~cellfun('isempty',aux)); % {'Rpu','Lpu'}

%=== Estimation of TargetFlow in systemic ArtVen's for flow control
qSy   = Get('ArtVen','q',Sy);% simulated flows in systemic ArtVen
qPu   = Get('ArtVen','q',Pu);% simulated flows in pulmonary ArtVen
q0AvSy= Get('ArtVen','q0AV',Sy); % Get reference flow per systemic ArtVen

if P.General.SaturationControl
    % 0: flow set / 1: control of flow by AV-saturation difference
    % Preparation execution of Saturation.m
    % Setting venous target saturation allowing execution of Saturation.m
    % Quantitative specific data in program text ++++++++++
    SatVenSy=[0.50,0.75,0.75,0.75]; %{'Ca','Br','Ce','Fe'}
    SatVenPu=0.98;
    Put('ArtVen','SatVen',Sy,SatVenSy);
    Put('ArtVen','SatVen',Pu,SatVenPu);
    P.General.nBeatSat=5; % number of cycles for to approach steady state of saturation
    Saturation
    
    % TargetFlow qSyTarget determined by systemic oxygen consumption
    % expressed as qDSatSyTarget= flow*DSaturation(A-V) in systemic
    % ArtVen's.
    % Note that DSaturation(A-V) ~= SatNodeProx - SatNodeDist
    SatSyProx      = P.Node.Sat(:,Get('ArtVen','iNode',Sy)); % proximal node saturation
    % qDSatSyTarget  = q0AvSy.*(mean(SatSyProx)-SatVenSy);
    % method to calculate qDSatSyTarget
    qDSatSyTarget  = [0.0722,0.0277,0.1316,0.0369]*P.General.q0;
    qDSatSy        = mean(bsxfun(@minus,SatSyProx,SatVenSy).*qSy); % mean flow*DavSaturation product
    FacqDSatControl= exp((1-qDSatSy./qDSatSyTarget)); % Vasodilation factor
    qSyTarget      = q0AvSy.*FacqDSatControl; % Target flow in systemic ArtVen's
    % Vasodilation if consumption > delivery by actual flow
    Put('ArtVen','q0AV',Sy,q0AvSy.*(qSyTarget./q0AvSy).^FbFac);
else % direct flow control, ratio of flows in systemic ArtVens preserved
    qSyTarget=q0AvSy/sum(q0AvSy)*P.General.q0; % Target flow in systemic ArtVen's
    Put('ArtVen','q0AV',Sy,qSyTarget);
end
% save Sy-Target flows in P.ArtVen.q0AV
%=== End of Estimation of TargetFlow

%Control of pressure by change of circulating volume
FacpControl= mean(P.Node.p(:,P.Node.iBaro))/P.General.p0;
P.General.FacpControl=FacpControl^FbFac; % >1 -> Volume decrease

% Control of flow to P.ArtVen.q0AV for systemic ArtVen's
p0AvSy=Get('ArtVen','p0AV',Sy);
kAvSy =Get('ArtVen','kAV' ,Sy);
FacqControlSy = exp(FbFac*(1-mean(qSy)./q0AvSy));
Put('ArtVen','p0AV',Sy,p0AvSy.*FacqControlSy.^(-1./kAvSy)...
    /P.General.FacpControl);
% Division by FacpControl improves stability of feedback

% Flow in Sys and Pu to judge steady state
FlowVec=[sum(mean(qSy)),sum(mean(qPu))]; % Systemic/Pu flow
disp('Flow/q0: Sys Pu');
disp(FlowVec/P.General.q0);

% Estimate AV-delay
P.General.TauAv=0.185*P.General.tCycle;

% ====== AdaptSpecific, different for Rest, Exercise
NoAdapt= strcmp(P.General.AdaptFunction,'Adapt0P');
if NoAdapt % No adaptation, regular sequence of beats
    % === Faster Steady State at rest
    VecV=Get('Cavity','V','All');
    P.General.In =[P.General.In ;VecV(  1,:) ];
    P.General.Out=[P.General.Out;VecV(end,:) ];
    if P.General.Fast;% Put new start values in structure P
        Vec= SteadyStateP;
        i=1  ; j=P.Cavity.n; Put('Cavity','V','All',Vec(i:j));
    end
else
    feval(P.General.AdaptFunction); % specific adapt function
end

% Judging quality of steady state
if size(P.General.Out,1)>1; % Escape if steady state is reached
    ErrVec= 1000*log( P.General.Out(end,:)./P.General.In(end,:) );
    disp(['Stationarity error= ',...
        num2str( round(sqrt(mean(ErrVec .^ 2 ))) )] );
    %=== ERROR criterium on flow stationarity
    if sqrt(mean(ErrVec .^ 2 ))<1 && (~NoAdapt || P.General.Fast)
        P.General.tEnd= P.General.tEnd-0.5*(P.General.tEnd-P.t(end));
    end
end;

% get the initial condition for next beat, P is most compact information to
% start the simulation
P2SVar; % load physiologic data in record of state variables P.SVar

end

