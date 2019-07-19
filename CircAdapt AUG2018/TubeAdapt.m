function TubeAdapt(StrTube,AdaptType)
% function TubeAdapt(StrTube,AdaptType);
% Adaptation of Diameter and Wall thickness of Tube to
% pressure and flow.
% StrTube= array of Tube names, e.g. {'AoBr','BrCe'}
% AdaptType= {'Diameter', 'WallVolume'} indicates type of adaptation
% Theo Arts, Maastricht University, April 2, 2014

global P
if P.Tube.n==0; return; end;

% Determine ArtVen indices
iTube= Get('Tube','Index',StrTube); %get related tube indices
Check=@(Str) ~isempty(strmatch(Str,AdaptType,'exact'));
if Check('All')
    AdaptDiameter=1; AdaptWallVolume=1;
else
    AdaptDiameter   =Check('Diameter');
    AdaptWallVolume =Check('WallVolume');
end

q     = (P.Tube.qProx(:,iTube)+P.Tube.qDist(:,iTube))/2;
A     = P.Tube.A(:,iTube);
pTrans= P.Tube.pTrans(:,iTube);
k     = P.Tube.k(:,iTube);
Z     = (P.Tube.ZL(:,iTube)+P.Tube.ZR(:,iTube))/2;

AWall = P.Tube.AWall(iTube);
A0    = P.Tube.A0(iTube);
p0    = P.Tube.p0(iTube);
vImpact   = P.Tube.Adapt.vImpact(iTube); % assumed body impact velocity
vFlowMean = P.Tube.Adapt.vFlowMean(iTube); % target value flow velocity
vMean = mean(abs(q./A)); %changed O402
vMean = max(mean(max(q,0)./A),mean(max(-q,0)./A)); %changed O402
pMean = sqrt(mean(pTrans).^2);

pMax= sqrt(1.7*mean(pTrans.^2)) + max(Z .* A) .* vImpact; % maximum pressure
AMax= max(A);
WallStress= 3*pMax.*(0.5+AMax./AWall); % maximum wall stress

% Correction of wall shear rate for small blood vessels
ASmall   = 0.7698e-6;% Diam=0.9 mm Ref Arts 2012
FacASmall= 1+sqrt(ASmall./A0);
Fac_vFlow    = FacASmall.*vMean./vFlowMean; % mean velocity/target value
FacWallStress= WallStress./P.Tube.Adapt.WallStress(:,iTube); % max wall stress/target value
Facp0        = pMean./p0; % mean pressure/target value
% FacA0        = AMean./A0; % mean cross-section/target value
%---

%=== Carrying out adaptation
FacV=ones(size(A0));
if AdaptDiameter; % adapts diameter to flow
    a= 0.20; % feedback factor, low value-> slow but stabile
    FacV = Fac_vFlow.^a; % volume stabilization during adaptation
    A0   = A0 .* Fac_vFlow.^a; % set working cross-section
    p0   = p0 .* Facp0.^(3*a./k); % set working pressure
end

if AdaptWallVolume % adapts wall thickness to pressure
    a=0.25; % feedback factor, low value-> slow but stabile
    AWall= AWall .* FacWallStress.^a;
end
%---
P.Tube.V(end,iTube)= P.Tube.V(end,iTube).*FacV; 
% SVar stabilization
P.Tube.A0(iTube)   =A0; % working cross-section
P.Tube.p0(iTube)   =p0; % working pressure
P.Tube.AWall(iTube)=AWall; % wall cross-section

aux=[P.Tube.Name;repmat({' '},[1,P.Tube.n])];
disp(['Deviation of ',aux{:},':'])
disp(['vFlow     : ', num2str(round(1000*log(Fac_vFlow    )),'%+6.3d')]);
disp(['WallStress: ', num2str(round(1000*log(FacWallStress)),'%+6.3d')]);

end

