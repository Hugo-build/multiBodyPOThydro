function plot_div_mode(pers, A_tot_WAMIT,B_tot_WAMIT, AB_tot, F_tot_WAMIT, F_tot, Ndl, path)

% selected added mass and damping plots (total)     
% AB_total
xlims = [0 5];
modestoplot = [1,1;2,2;3,3]; 
ff1 = figure;
set(ff1,'Position',[100 100 800 500])
for ii = 1:length(modestoplot(:,1))
    i1 = modestoplot(ii,1);
    i2 = modestoplot(ii,2);
    subplot(3,2,2*(ii-1)+1)
    plot(2*pi./pers,squeeze(A_tot_WAMIT(:,i1,i2)),'k','LineWidth',1.2)
    hold on
    for idl = 1:Ndl
        plot(2*pi./pers,squeeze(real(AB_tot(idl,:,i1,i2))),'--','LineWidth',1.2)
        hold on
    end
    ylabel(['A_' num2str(i1) '_' num2str(i2) ', kg' ])
    xlabel('\omega, rad/s')
    grid on
    legend('rigid body','dm1','dm2','dm3','dm4','dm5')
    xlim(xlims)
    
    subplot(3,2,2*ii)
    plot(2*pi./pers,squeeze(B_tot_WAMIT(:,i1,i2)),'k','LineWidth',1.2)
    hold on
    for idl = 1:Ndl
       plot(2*pi./pers,squeeze(2*pi./pers.*imag(AB_tot(idl,:,i1,i2))'),'--','LineWidth',1.2)
       hold on
    end
    ylabel(['B_' num2str(i1) '_' num2str(i2) ', kg'  '/s'])
    xlabel('\omega, rad/s')
    grid on
    legend('rigid body','dm1','dm2','dm3','dm4','dm5')
    xlim(xlims)
end
saveas(ff1,[path 'AB_ii_check.png'])
%--------------------------------------------------------------------------
ff2 = figure;
set(ff2,'Position',[100 100 800 500])
modestoplot = [1,5;2,4;2,6];
for ii = 1:length(modestoplot(:,1))
    i1 = modestoplot(ii,1);
    i2 = modestoplot(ii,2);
    subplot(3,2,2*(ii-1)+1)
    plot(2*pi./pers,squeeze(A_tot_WAMIT(:,i1,i2)),'k','LineWidth',1.2)
    hold on
    for idl = 1:Ndl
        plot(2*pi./pers,squeeze(real(AB_tot(idl,:,i1,i2))),'--','LineWidth',1)
        hold on
    end
    ylabel(['A_' num2str(i1) '_' num2str(i2) ', kgm'  ])
    xlabel('\omega, rad/s')
    grid on
    legend('rigid body','dm1','dm2','dm3','dm4','dm5')
    xlim(xlims)
    
    subplot(3,2,2*ii)
    plot(2*pi./pers,squeeze(B_tot_WAMIT(:,i1,i2)),'k','LineWidth',1.2)
    hold on
    for idl = 1:Ndl
       plot(2*pi./pers,squeeze(2*pi./pers.*imag(AB_tot(idl,:,i1,i2))'),'--','LineWidth',1)
       hold on
    end
    ylabel(['B_' num2str(i1) '_' num2str(i2) ', kgm' '/s'])
    xlabel('\omega, rad/s')
    grid on
    legend('rigid body','dm1','dm2','dm3','dm4','dm5')
    xlim(xlims)
end
saveas(ff2,[path 'AB_ij_check.png'])
%--------------------------------------------------------------------------
ff3 = figure;
set(ff3,'Position',[100 100 800 500])
modestoplot = [5,5;4,4;6,6];
for ii = 1:length(modestoplot(:,1))
    i1 = modestoplot(ii,1);
    i2 = modestoplot(ii,2);
    subplot(3,2,2*(ii-1)+1)
    plot(2*pi./pers,squeeze(A_tot_WAMIT(:,i1,i2)),'k','LineWidth',1.2)
    hold on
    for idl = 1:Ndl
        plot(2*pi./pers,squeeze(real(AB_tot(idl,:,i1,i2))),'--','LineWidth',1)
        hold on
    end
    ylabel(['A_' num2str(i1) '_' num2str(i2) ', kgm^2'  ])
    xlabel('\omega, rad/s')
    grid on
    legend('rigid body','dm1','dm2','dm3','dm4','dm5')
    xlim(xlims)
    
    subplot(3,2,2*ii)
    plot(2*pi./pers,squeeze(B_tot_WAMIT(:,i1,i2)),'k','LineWidth',1.2)
    hold on
    for idl = 1:Ndl
       plot(2*pi./pers,squeeze(2*pi./pers.*imag(AB_tot(idl,:,i1,i2))'),'--','LineWidth',1)
       hold on
    end
    ylabel(['B_' num2str(i1) '_' num2str(i2) ',kgm^2 ' '/s'])
    xlabel('\omega, rad/s')
    grid on
    legend('rigid body','dm1','dm2','dm3','dm4','dm5')
    xlim(xlims)
end
saveas(ff3,[path 'AB_jj_check.png'])
%--------------------------------------------------------------------------
ff4 = figure;
set(ff4,'Position',[100 100 800 500])
modestoplot = [1,2;3,4;5,6];
for ii = 1:length(modestoplot(:,1))
    i1 = modestoplot(ii,1);
    i2 = modestoplot(ii,2);
    subplot(3,4,4*(ii-1)+1)
    plot(2*pi./pers,squeeze(abs(F_tot_WAMIT(:,1,i1))'),'k','LineWidth',2)
    hold on
    for idl = 1:Ndl
        plot(2*pi./pers,squeeze(abs(F_tot(idl,:,1,i1))'),'LineWidth',2)
    end
    
    subplot(3,4,4*(ii-1)+2)
    plot(2*pi./pers,2*pi/180*squeeze(angle(F_tot_WAMIT(:,1,i1))'),'k','LineWidth',2)
    hold on
    for idl = 1:Ndl
        scatter(2*pi./pers,2*pi/180*squeeze(angle(F_tot(idl,:,1,i1))'),'LineWidth',2)
    end
    
    subplot(3,4,4*(ii-1)+3)
    plot(2*pi./pers,squeeze(abs(F_tot_WAMIT(:,8,i2))'),'k','LineWidth',2)
    hold on
    for idl = 1:Ndl
        plot(2*pi./pers,squeeze(abs(F_tot(idl,:,8,i2))'),'LineWidth',2)
    end
    
    subplot(3,4,4*ii)
    plot(2*pi./pers,2*pi/180*squeeze(angle(F_tot_WAMIT(:,8,i2))'),'k','LineWidth',2)
    hold on
    for idl = 1:Ndl
        scatter(2*pi./pers,2*pi/180*squeeze(angle(F_tot(idl,:,8,i2))'),'LineWidth',2)
    end

end
saveas(ff4,[path 'F_jj_check.png'])






