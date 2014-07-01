function plotseg(seg,raw,interpo,smooth,bary)
    colors = {'.y','.g','.b','.m','.c','.k'};
    for i = 1:length(seg)
        if interpo
            plot(seg{i}.xq, seg{i}.yq, colors{rem(i+1,6)+1}), hold on
        end
        if raw
            plot(seg{i}.xi,seg{i}.yi,colors{rem(i,6)+1}), hold on
        end
        if smooth
            plot(seg{i}.xs, seg{i}.ys,'-r','LineWidth',2), hold on
        end
        if bary
            plot(seg{i}.xb, seg{i}.yb, '.k','MarkerSize',20), hold on
        end
    end
    xlim([0 78]),ylim([0 101]), hold off
end