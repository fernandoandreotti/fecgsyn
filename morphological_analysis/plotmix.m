%% Function to plot generated mixtures
function plotmix(out)
FONT_SIZE = 15;
LINE_WIDTH = 1.2;
MARKER_SIZE = 5;
N = 5;      % number of plots
Nchan = size(out.mixture,1);
chan = randsample(1:Nchan,N);
% watch out! tight_subplot has a different license and therefore should not
% be open-sourced with our code.
figure(1)
clf
h = tight_subplot(N,1,[.01 .03],[.05 .05],[.05 .05]);

for i = 1:N
    axes(h(i))
    plot(out.mixture(chan(i),:),'Color',[.6 .6 .6],'LineWidth',LINE_WIDTH)
    hold on
    for j = 1:length(out.fqrs)
        valf = 1.2*median(out.mixture(chan(j),out.fqrs{j}))*ones(1,length(out.fqrs{j}));
        plot(out.fqrs{j},valf,'xr','MarkerSize',MARKER_SIZE,'MarkerFaceColor','r')
    end
    valm = 1.2*median(out.mixture(chan(i),out.mqrs))*ones(1,length(out.mqrs));
    plot(out.mqrs,valm,'dk','MarkerSize',MARKER_SIZE,'MarkerFaceColor','k')   
    ylabel(['ch' num2str(chan(i))],'FontSize',FONT_SIZE)
    hold off
end
linkaxes(h)
