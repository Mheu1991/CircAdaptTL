figure(2); plot(sum(P.Cavity.V,2)+sum(P.Tube.V,2));
V=sum(mean(P.Cavity.VDot));
T=sum(mean(P.Tube.VDot));
disp('Vol Cavities, Tubes, Total (ml)')
disp(1e6*[V,T,V+T])

