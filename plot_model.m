function plot_model(idl,crd_C, crd_Cs, zC, lengthP, crd_Pon, crd_P1, crd_P2, crd_P3, xP1, XCT, YCT, ZCT, nodes,...
                    Secs,SecP1, SecP2, SecP3, SecC, SecCs1, SecCs2, SecCs3)
%% -------------------------------------------
% plot cross sections
close all
path = 'multibody_total\';
%---------------------------------------------------
% column sections
figure
subplot(1,2,1)
scatter3(crd_C(:,1),crd_C(:,2),crd_C(:,3), 50 ,[0.7 0.7 0.7],'square')
axis([-40 60 -60 60 -20 15])
hold on
scatter3(crd_Cs(:,1),crd_Cs(:,2),crd_Cs(:,3), 50 ,[0.7 0.7 0.7],'square')
hold on
x = -8:1:8; y = -8:1:8; z = -20:1:15; [x,y,z] = meshgrid(x,y,z); v = x*0 +1;
xslice = []; yslice = []; zslice = zC(2:end-1)';
slice(x,y,z,v,xslice,yslice,zslice)
hold on
x = x + lengthP;
slice(x,y,z,v,xslice,yslice,zslice)
hold on
x = x - lengthP - lengthP*cosd(60); y1 = y+lengthP *sind(60); y2 = y+lengthP *sind(-60);
slice(x,y1,z,v,xslice,yslice,zslice)
hold on
slice(x,y2,z,v,xslice,yslice,zslice)

%--------------------------------------------------------
% pontoon sections
subplot(1,2,2)
scatter3(crd_Pon(:,1),crd_Pon(:,2),crd_Pon(:,3), 50 ,[0.7 0.7 0.7],'square')
axis([-40 60 -60 60 -20 0])
hold on
x = -40:1:60; y = -8:1:8; z = -20:1:-13;
[x,y,z] = meshgrid(x,y,z); v = x*0 +1;
xslice = xP1(1:end-1)'; yslice = []; zslice = [];
s1 = slice(x,y,z,v,xslice,yslice,zslice); s2 = slice(x,y,z,v,xslice,yslice,zslice); s3= slice(x,y,z,v,xslice,yslice,zslice);
rotate(s2, [0 0 1], 120, [0, 0, 0]); rotate(s3, [0 0 1], -120, [0, 0, 0]);

%----------------------------------------------------------------------------------------------------
% Riflex beam
ff1 = figure;
scatter3(XCT, YCT, ZCT, 50 ,[0.9 0.9 0.9],'square')
hold on
scatter3(Secs(:,1),Secs(:,2),Secs(:,3),...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r')
axis([-40 60 -60 60 -20 15])
hold on; plot3(SecP1(:,1),SecP1(:,2),SecP1(:,3),'color',[0.4660 0.6740 0.1880],'linewidth',4)
hold on; plot3(SecP2(:,1),SecP2(:,2),SecP2(:,3),'color',[0.4660 0.6740 0.1880],'linewidth',4)
hold on; plot3(SecP3(:,1),SecP3(:,2),SecP3(:,3),'color',[0.4660 0.6740 0.1880],'linewidth',4)
hold on; plot3(SecCs1(:,1),SecCs1(:,2),SecCs1(:,3),'color',[0.4660 0.6740 0.1880],'linewidth',4)
hold on; plot3(SecCs2(:,1),SecCs2(:,2),SecCs2(:,3),'color',[0.4660 0.6740 0.1880],'linewidth',4)
hold on; plot3(SecCs3(:,1),SecCs3(:,2),SecCs3(:,3),'color',[0.4660 0.6740 0.1880],'linewidth',4)
hold on; plot3(SecC(:,1),SecC(:,2),SecC(:,3),'color',[0.4660 0.6740 0.1880],'linewidth',4)
saveas(ff1,[path 'R_beam' num2str(idl) '.png'])

%----------------------------------------------------------------------------------------------------
% SIMO body
ff2 = figure;
scatter3(XCT, YCT, ZCT, 50 ,[0.53 0.81 0.92],'square')
hold on
scatter3(nodes(:,1),nodes(:,2),nodes(:,3),...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','b')
axis([-40 60 -60 60 -20 15])
%---------------------------------------------------------------------------------------------------

figure
scatter3(crd_P1(:,1),crd_P1(:,2),crd_P1(:,3))
axis([-40 60 -60 60 -20 0])
xlabel('X [m]')
ylabel('Y [m]')
hold on
scatter3(crd_P2(:,1),crd_P2(:,2),crd_P2(:,3),'r')
hold on
scatter3(crd_P3(:,1),crd_P3(:,2),crd_P3(:,3),'k')
saveas(ff2,[path 'S_body' num2str(idl) '.png'])

%--------------------------------------------------------------------------------------------------------------------
% %%
% figure
% scatter3(XCT, YCT, ZCT)
% 
% figure
% crd_Pon = [crd_P1;crd_P2;crd_P3];
% scatter3(crd_Pon(:,1),crd_Pon(:,2),crd_Pon(:,3))
% axis([-40 60 -60 60 -20 0])
% 
% figure
% % scatter3(crd_BN(:,1),crd_BN(:,2),crd_BN(:,3))
% % axis([-40 60 -60 60 -20 0])
% scatter(crd_BN(:,1),crd_BN(:,2))
% axis([-40 60 -60 60])
% 
% figure
% scatter3(crd_C(:,1),crd_C(:,2),crd_C(:,3))
% axis([-40 60 -60 60 -20 0])
% 
% figure
% crd_Cs = [crd_Cs1;crd_Cs2;crd_Cs3];
% scatter3(crd_Cs(:,1),crd_Cs(:,2),crd_Cs(:,3))
% axis([-40 60 -60 60 -20 0])
