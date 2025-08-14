function plot_excitation(pers, F_tot_WAMIT, F_tot)


% excitation     % F_tot( period, direction, 6)
xlims = [0 3];
ff = figure;
set(ff,'Position',[100 100 800 500])
for ii = [1,3,5]
    subplot(3,2,ii)
    plot(2*pi./pers,squeeze(abs(F_tot_WAMIT(:,1,ii))),'k')
    hold on
    plot(2*pi./pers,squeeze(abs(F_tot(:,1,ii))),'r--')
    
    legend('cyl - rigid body','cyl - distributed')
   
    unitname = 'N/m';
    if ii>=4; unitname = 'Nm/m'; end
    ylabel(['|X_' num2str(ii) '|, ' unitname ])
    grid on
    xlabel('\omega, rad/s')
    xlim(xlims)
      
    subplot(3,2,ii+1)
    plot(2*pi./pers,squeeze(abs(F_tot_WAMIT(:,8,ii+1))),'k')
    hold on
    plot(2*pi./pers,squeeze(abs(F_tot(:,8,ii+1))),'r--')
    ylabel(['phase(X_' num2str(ii) '), deg' ])
    grid on
    xlabel('\omega, rad/s')
   
    legend('cyl - rigid body','cyl - distributed')
    unitname = 'N/m';
    if ii>3; unitname = 'Nm/m'; end
    ylabel(['|X_' num2str(ii) '|, ' unitname ])
    grid on
    xlabel('\omega, rad/s')
    xlim(xlims)
end