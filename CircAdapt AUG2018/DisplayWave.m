% Script DisplayWave
% Theo Arts, Maastricht University, Oct 13, 2012
% Displays wave pressures and flow velocities, based on structure P

pAr={'Lv','Ao','BrAr','CeAr','FeAr'};
pVe={'Ra','Vc','BrVe','CeVe','FeVe'};
TubeArtFlow={'AoBrAr','BrArCeAr','CeArFeAr'};
TubeVenFlow={'VcBrVe','BrVeCeVe','CeVeFeVe'};
ArtFlowDist={'CaAr','BrAr','CeAr','FeAr'};

pThorax=P.Tube.p(:,1)-P.Tube.pTrans(:,1);
pPeri  =Get('Bag','p','Peri');
T=P.t-P.t(1);

figure(2); M=Get('Node','p',pAr);
plot(T,M);
axis([min(T),max(T),min(M(:)),max(M(:))])
legend(pAr)
title('Arterial Pressures')

figure(3); M=[Get('Node','p',pVe),pThorax,pPeri];
plot(T,M);
axis([min(T),max(T),min(M(:)),max(M(:))])
legend([pVe,{'Thorax','Peri'}])
title('Venous Pressures')

M1= [Get('Tube','qProx',TubeArtFlow),Get('Tube','qDist',TubeArtFlow)];
A= Get('Tube','A',TubeArtFlow);
M=M1./[A,A];
figure(4);
plot(T,M);
axis([min(T),max(T),min(M(:)),max(M(:))])
legend([TubeArtFlow,TubeArtFlow])
title('Prox and Dist Artery Velocity')

M1= [-Get('Tube','qProx',TubeVenFlow),-Get('Tube','qDist',TubeVenFlow)];
A= Get('Tube','A',TubeVenFlow);
M= M1./[A,A];
figure(5);
plot(T,M);
axis([min(T),max(T),min(M(:)),max(M(:))])
legend([TubeVenFlow,TubeVenFlow])
title('Prox and Dist Venous Velocity')

iCav=Get('ArtVen','iCavity','All');
Len= Get('ArtVen','Len','All');
nt=length(P.t);
AAr=P.Cavity.V(:,iCav)./repmat(Len(1,:),[nt,1]);
M1=Get('ArtVen','q','All');
M2=P.Cavity.VDot(:,Get('ArtVen','iCavity','All'));
M=(M1+M2)./AAr;
figure(6);
plot(T,M);
axis([min(T),max(T),min(M(:)),max(M(:))])
legend(P.ArtVen.Name)
title('Arterial Flow Velocity into Periphery')

