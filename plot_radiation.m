function plot_radiation(pers, A_tot_WAMIT,B_tot_WAMIT, AB_tot)

% selected added mass and damping plots (total)     
% AB_total
xlims = [0 3];
modestoplot = [1,1;2,2;3,3]; 
unitnames = {'kg','kg','kgm^2','kgm'};
ff = figure;
%set(ff,'Position',[100 100 900 1200])
for ii = 1:length(modestoplot(:,1))
    i1 = modestoplot(ii,1);
    i2 = modestoplot(ii,2);
    subplot(3,2,2*(ii-1)+1)
    plot(2*pi./pers,squeeze(A_tot_WAMIT(:,i1,i2)),'k')
    hold on
    plot(2*pi./pers,squeeze(real(AB_tot(:,i1,i2))),'r--')
    ylabel(['A_' num2str(i1) '_' num2str(i2) ', ' unitnames{ii} ])
    xlabel('\omega, rad/s')
    grid on
    if ii==1; legend('cyl - rigid body','cyl - distributed'); end
    xlim(xlims)
    
    subplot(3,2,2*ii)
    plot(2*pi./pers,squeeze(B_tot_WAMIT(:,i1,i2)),'k')
    hold on
    plot(2*pi./pers,squeeze(-2*pi./pers.*imag(AB_tot(:,i1,i2))),'r--')
    ylabel(['B_' num2str(i1) '_' num2str(i2) ', ' unitnames{ii} '/s'])
    xlabel('\omega, rad/s')
    grid on
    xlim(xlims)
end
%--------------------------------------------------------------------------
figure
modestoplot = [1,5;2,4;2,6];
for ii = 1:length(modestoplot(:,1))
    i1 = modestoplot(ii,1);
    i2 = modestoplot(ii,2);
    subplot(3,2,2*(ii-1)+1)
    plot(2*pi./pers,squeeze(A_tot_WAMIT(:,i1,i2)),'k')
    hold on
    plot(2*pi./pers,squeeze(real(AB_tot(:,i1,i2))),'r--')
    ylabel(['A_' num2str(i1) '_' num2str(i2) ', ' unitnames{ii} ])
    xlabel('\omega, rad/s')
    grid on
    if ii==1; legend('cyl - rigid body','cyl - distributed'); end
    xlim(xlims)
    
    subplot(3,2,2*ii)
    plot(2*pi./pers,squeeze(B_tot_WAMIT(:,i1,i2)),'k')
    hold on
    plot(2*pi./pers,squeeze(-2*pi./pers.*imag(AB_tot(:,i1,i2))),'r--')
    ylabel(['B_' num2str(i1) '_' num2str(i2) ', ' unitnames{ii} '/s'])
    xlabel('\omega, rad/s')
    grid on
    xlim(xlims)
end
%--------------------------------------------------------------------------
figure
modestoplot = [5,5;4,4;6,6];
for ii = 1:length(modestoplot(:,1))
    i1 = modestoplot(ii,1);
    i2 = modestoplot(ii,2);
    subplot(3,2,2*(ii-1)+1)
    plot(2*pi./pers,squeeze(A_tot_WAMIT(:,i1,i2)),'k')
    hold on
    plot(2*pi./pers,squeeze(real(AB_tot(:,i1,i2))),'r--')
    ylabel(['A_' num2str(i1) '_' num2str(i2) ', ' unitnames{ii} ])
    xlabel('\omega, rad/s')
    grid on
    if ii==1; legend('cyl - rigid body','cyl - distributed'); end
    xlim(xlims)
    
    subplot(3,2,2*ii)
    plot(2*pi./pers,squeeze(B_tot_WAMIT(:,i1,i2)),'k')
    hold on
    plot(2*pi./pers,squeeze(-2*pi./pers.*imag(AB_tot(:,i1,i2))),'r--')
    ylabel(['B_' num2str(i1) '_' num2str(i2) ', ' unitnames{ii} '/s'])
    xlabel('\omega, rad/s')
    grid on
    xlim(xlims)
end