function plotspectrowmat(F,T,P)
    surf(T,F,10*log10(P),'edgecolor','none'); axis tight; view(0,90);
end