function ArtVenAdapt(StrAV,AdaptType)
% function ArtVenAdapt(StrAV,AdaptType);
% Adaptation of Diameter and Wall thickness of Art and Ven to
% pressure and flow.
% StrAV= array of ArtVen names, e.g. {'Br','Pu'}
% AdaptType= {'Diameter', 'WallVolume'} indicates type of adaptation
% Theo Arts, Maastricht University, Oct 13, 2012

global P

FeedBack= 0.15; % feedback factor, low value-> slow but stabile

% Determine ArtVen indices
iAV   = Get('ArtVen','Index',StrAV); % read string -> index
nAV   = size(iAV,2);

% read type(s) of adaptation
iAdapt= Str2Index(AdaptType,{'Diameter','WallVolume'});
Adapt=[0,0]; Adapt(iAdapt)=1;
AdaptDiameter=Adapt(1); AdaptWallVolume=Adapt(2); % active if Adapt==1

szAV=[2,nAV]; szRow=[1,2*nAV]; % shape of array in P.ArtVen and P.Cavity
iCavArt = P.ArtVen.iCavity(:,iAV); iCavVen=iCavArt+1; 
iCav    = reshape([iCavArt;iCavVen],szRow); % corresponding Cavity indices in row
q     = P.ArtVen.q(:,reshape([iAV;iAV],szRow)); % Microcirculatory flow
A     = P.Cavity.A(:,iCav); % cavity cross-section
pTrans= P.Cavity.pTrans(:,iCav); % cavity pressure
ArtVenV2p;% recalculates wave impedance before waterfall corrections
Z     = P.Cavity.Z(:,iCav); % cavity wave impedance

vMean= reshape(max(mean(max(q,0)./A),mean(max(-q,0)./A)),szAV);
% mean flow velocity O402

pMean= reshape(sqrt( mean(pTrans.^2)),szAV); % mean transmural pressure

AWall  = P.ArtVen.AWall(:,iAV); % wall cross-section
p0     = P.ArtVen.p0(:,iAV); % working pressure
A0     = P.ArtVen.A0(:,iAV); % working cross-section
k      = P.ArtVen.k(:,iAV); % working cross-section
vImpact   = P.ArtVen.Adapt.vImpact(:,iAV); % assumed body impact velocity
vFlowMean = P.ArtVen.Adapt.vFlowMean(:,iAV); % target value flow velocity

pMax = max(pTrans) + max(bsxfun(@times,Z.*A, vImpact(:)')); % maximum pressure
pMax = reshape(pMax,szAV);
AMax = max(A); % maximum cross-sectional area
AMax = reshape(AMax,szAV);
WallStress= 3*pMax.*(0.5+AMax./AWall); % maximum wall stress

% For small blood vessels shear rate determines diameter (1/r3), for large
% ones mean blood velocity is controlled (1/r2)
ASmall   = 0.7698e-6;% Critical diameter=0.9 mm Ref Arts 2012
FacASmall= 1+sqrt(ASmall./A0); % for small A0 constant wall shear rate
Fac_vFlow    = FacASmall.*vMean./vFlowMean; % mean velocity/target value
FacWallStress= WallStress./P.ArtVen.Adapt.WallStress(:,iAV); % max wall stress/target value
Facp0        = max(pMean,0.1*pMax)./p0; % mean pressure/target value
% For extremely low pMean, a minimum pressure level is used =0.1 max(p)
%---

%=== Carrying out adaptation
FacV=ones(size(A0)); % blood volume correction after geometric adaptation
if AdaptDiameter; % adapts diameter to flow
    FacV = reshape(Fac_vFlow.^FeedBack,szRow); % volume stabilization during adaptation
    A0 = A0 .* Fac_vFlow.^FeedBack; % cross-section adapts to fixed mean velocity
    p0 = p0 .* Facp0.^(3*FeedBack./k); % p0 set near pMean, but not equal
end
if AdaptWallVolume % adapts wall thickness to pressure
    AWall= AWall .* FacWallStress.^FeedBack;
end
%---

P.Cavity.V(end,iCav)= P.Cavity.V(end,iCav).*FacV(:)'; 
% SVar-volumes adapted so that pressure is maintained
P.ArtVen.A0(:,iAV)   =A0; % working cross-section
P.ArtVen.p0(:,iAV)   =p0; % working pressure
P.ArtVen.AWall(:,iAV)=AWall; % wall cross-section

aux=[P.ArtVen.Name;repmat({' '},[1,P.ArtVen.n])];
disp(['Deviation of ',aux{:},':'])
disp(['vFlow     : ', num2str(round(1000*log(Fac_vFlow(:)'    )),'%+6.3d')]);
disp(['WallStress: ', num2str(round(1000*log(FacWallStress(:)')),'%+6.3d')]);

end

function Index=Str2Index(Str,Names)
% conversion of Name(s) to indices in structure
n=size(Names,2); iN=(1:n)';
Index=[];
if ischar(Str) %if single string
    if strcmp(Str,'All')
        Index=1:n;
    else
        Index=strcmp(Str,Names)*iN;
    end
else % if array of strings
    for i=1:length(Str)
        Index= [Index,strcmp(Str{i},Names)*iN];
    end
end
end

