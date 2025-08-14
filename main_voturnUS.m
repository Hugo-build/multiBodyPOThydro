%% Create Wamit results for sections
% Created by Xiaoming. reference including Spar code and Svendsen 2016
% 01.10.2021
% only translational dof (both calculation and output), for SIMO body way


% contents
% 1-input,
% 2-Split TLP and define multibody, (section number)
% 3-Generate SIMO body, integration
% 4-write results to file. 

%             dl = 1	 dl = 2	   dl = 3	 dl = 4   dl =5
% Column x4	    16	        12	      8	       4        2
% Pontoon x3	16	        12	      8	       4        2
% Total	        108	        84	      56	   28       14
%% Section 1-Input %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
tic
set(0,'defaultfigurecolor','w')
format long
%% output and input
prefix = 'semi';              
folder = 'Wam_semi\';      
folderout='matlab_out_4April4_2\';                
mkdir matlab_out_4April4_2             % output
mode_out = 3;                    % 3 for potential flow library, 6 for SIMO body way. model for output

%% Read WAMIT output
datGDF = dlmread([folder prefix '.GDF'],'',1,0);    % geometric data file
datpnl = dlmread([folder prefix '.pnl']);           % panel data file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
dat1 = dlmread([folder prefix '.1']);               % added mass and damping coefficients
dat2 = dlmread([folder prefix '.2']);               % exciting forces from Haskind relations
dat3 = dlmread([folder prefix '.3']);               % exciting forces from diffraction potential
dat5p = dlmread([folder prefix '.5p']);             % Hydrodynamic pressure on body surface

%% Define basic parameters
ULEN = datGDF(1,1);                   
g = datGDF(1,2);                      
XV = datGDF(4:end,1);                
YV = datGDF(4:end,2);
ZV = datGDF(4:end,3);

Mpnl = datpnl(:,1);                    
Kpnl = datpnl(:,2);                  
XCT = datpnl(:,3);  % x coordinate of panels               
YCT = datpnl(:,4);  % y coordinate of panels    
ZCT = datpnl(:,5);  % z coordinate of panels   
AREA = datpnl(:,6);
nx = datpnl(:,7);                      
ny = datpnl(:,8); 
nz = datpnl(:,9);
nrx = datpnl(:,10);                    
nry = datpnl(:,11);
nrz = datpnl(:,12);


pers = flipud(unique(dat5p(:,1)));             
Npers = length(pers);                          % number of periods
Npanels = length(Kpnl);                        % number of panels
Nbetas = length(unique(dat3(:,2)));            % number of wave headings
betas = unique(dat5p(Npanels+1:Npanels+Npanels*Nbetas,2));  
omega = 2*pi./pers;        
rho = 1025;
L = ULEN;                  

x0=0; %Position of new local coordinate system   
y0=0; %Position of new local coordinate system


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% section 2.   %%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Split Semi into column, base node (lower column) and (indivdual) pontoons
heightC = 13.0;                % upper column length (under water) 
heightPon = 7;
lengthP = 51.75; 

radC_c = 10/2;                   % central column diameter
radC_s = 12.5/2;                 % side column diameter
draft = min(ZCT);                % semi-sub draft
lengthBN = -draft-heightC;       % Base Node means BN
Zs = draft+8.5*0.5;              % spokes centerline in z, shouled be -

crd_Semi = [XCT YCT ZCT];            % all panel center coordinates Semi

pnl_C = find( sqrt(crd_Semi(:,1).^2 + crd_Semi(:,2).^2) <= radC_c & crd_Semi(:,3)>-heightC ); 
crd_C = crd_Semi(pnl_C,:);           % panel coordinates central column, x y z all

pnl_Cs1 = find( sqrt(crd_Semi(:,1).^2 + crd_Semi(:,2).^2) > radC_c & crd_Semi(:,3) > -heightC & crd_Semi(:,1) < 0 ); 
crd_Cs1 = crd_Semi(pnl_Cs1 ,:);      % panel coordinates side column 1
pnl_Cs2 = find( sqrt(crd_Semi(:,1).^2 + crd_Semi(:,2).^2) > radC_c & crd_Semi(:,3) > -heightC & crd_Semi(:,2) > 0 & crd_Semi(:,1) > 0 );
crd_Cs2 = crd_Semi(pnl_Cs2 ,:);      % panel coordinates side column 2
pnl_Cs3 = find( sqrt(crd_Semi(:,1).^2 + crd_Semi(:,2).^2) > radC_c & crd_Semi(:,3) > -heightC & crd_Semi(:,2) < 0 & crd_Semi(:,1) > 0 );
crd_Cs3 = crd_Semi(pnl_Cs3 ,:);      % panel coordinates side column 3

pnl_P1 = find(  crd_Semi(:,3) <= -heightC & crd_Semi(:,1) < -4 );   % x positive, 
crd_P1 = crd_Semi(pnl_P1,:);         % panel coordinates pontoon 1
pnl_P2 = find(  crd_Semi(:,3) <= -heightC & crd_Semi(:,2) > 0 & crd_Semi(:,1) > -4);
crd_P2 = crd_Semi(pnl_P2,:);         % panel coordinates pontoon 2
pnl_P3 = find(  crd_Semi(:,3) <= -heightC & crd_Semi(:,2) < 0 & crd_Semi(:,1) > -4);
crd_P3 = crd_Semi(pnl_P3,:);         % panel coordinates pontoon 3 % z+ r + xy 

% pnl_P1 = find( sqrt(crd_Semi(:,1).^2 + crd_Semi(:,2).^2) > radC_c & crd_Semi(:,3) <= -heightC & crd_Semi(:,1) > 4 );   % x positive, 
% crd_P1 = crd_Semi(pnl_P1,:);         % panel coordinates pontoon 1
% pnl_P2 = find( sqrt(crd_Semi(:,1).^2 + crd_Semi(:,2).^2) > radC_c & crd_Semi(:,3) <= -heightC & crd_Semi(:,2) > 0 & crd_Semi(:,1) < 4);
% crd_P2 = crd_Semi(pnl_P2,:);         % panel coordinates pontoon 2
% pnl_P3 = find( sqrt(crd_Semi(:,1).^2 + crd_Semi(:,2).^2) > radC_c & crd_Semi(:,3) <= -heightC & crd_Semi(:,2) < 0 & crd_Semi(:,1) < 4);
% crd_P3 = crd_Semi(pnl_P3,:);         % panel coordinates pontoon 3 % z+ r + xy 

pnl_BN = find( sqrt(crd_Semi(:,1).^2 + crd_Semi(:,2).^2) <= radC_c & crd_Semi(:,3) <= -heightC); 
crd_BN = crd_Semi(pnl_BN,:);           % panel coordinates central column, x y z all
% here, all panels are included. no missing, checked, correct
n_should =[length(pnl_C); length(pnl_Cs1); length(pnl_Cs2),; length(pnl_Cs3); length(pnl_P1); length(pnl_P2); length(pnl_P3)];
if Npanels ==sum(n_should)
   disp('all panels are included. no missing, checked, correct');
end
%% Define multi-body model
Ndl = 5;
idl = input('mode of sections divided == ');
if idl == 1
    NN = 16;
elseif idl == 2
    NN = 12;
elseif idl == 3
    NN = 8;
elseif idl == 4
    NN = 4;
elseif idl == 5
    NN = 2;
end
NsectionC = NN;    % number of section of Upper Column. 
NsectionBN = 1;   % number of section. 
NsectionP = NN;    % number of section. 

d1 = heightC/NsectionC;  %length of individual section        
d2 = lengthBN/NsectionBN;    
%--------------------------------------------------------------------------
% columns from low to top, the boundary of each section, not the center
zC = linspace(-heightC,0,NsectionC+1)';  
xC = zeros(size(zC)); % node x-locations on column
yC = zeros(size(zC)); % node y-locations on column
SecC = [xC yC zC];
% side columns  1 is along the x positive, 2 is in the 2nd qudarant. 3 is in 3rd.
zCs1 = zC ; zCs2 = zC ; zCs3 = zC ;
xCs1 = ones(size(zC))*lengthP;  xCs2 = ones(size(zC))*lengthP *(- cosd(60)); xCs3 = ones(size(zC))*lengthP *(- cosd(60));
yCs1 = zeros(size(zC));         yCs2 = ones(size(zC))*lengthP *(sind(60));   yCs3 = ones(size(zC))*lengthP *(- sind(60));
SecCs1 = [xCs1 yCs1 zCs1]; SecCs2 = [xCs2 yCs2 zCs2]; SecCs3 = [xCs3 yCs3 zCs3];
% pontoon 1 is along the x positive, 2 is in the 2nd qudarant. 3 is in 3rd.
xP1 = linspace(radC_c,radC_s+lengthP,NsectionP+1)';  % from column surface to tip. mesh number is 30 (29.7).
yP1 = zeros(size(xP1));                 % node y-locations on pontoon 1  1st 
zP1 = Zs*ones(size(xP1));               % node z-locations on pontoon 1
xP2 = -cosd(60)*xP1;                    % node x-locations on pontoon 2
yP2 = cosd(30)*xP1;                     % node y-locations on pontoon 2, 2nd qudarant
zP2 = Zs*ones(size(yP2));               % node z-locations on pontoon 2
xP3 = xP2;                              % node x-locations on pontoon 3   
yP3 = -yP2;                             % node y-locations on pontoon 3  3rd 
zP3 = zP2;                              % node z-locations on pontoon 3
SecP1 = [xP1 yP1 zP1]; SecP2 = [xP2 yP2 zP2]; SecP3 = [xP3 yP3 zP3];
Secs = [SecC; SecCs1; SecCs2; SecCs3; SecP1; SecP2; SecP3];
%---------objectives-----------------------
% SecC, SecCs1, SecCs2, SecCs3
% SecP1, SecP2, SecP3
% Secs
%--------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
% node, middle of section
% central column
zCn = zC - 0.5*d1*ones(size(zC));      
zCn = zCn(2:end);          
xCn = zeros(size(zCn));     
yCn = zeros(size(zCn));     
nodesC = [xCn yCn zCn];     
NnodesC = size(nodesC,1);  
% side column
zCsn1 = zCn ; zCsn2 = zCn ; zCsn3 = zCn ;
xCsn1 = ones(size(zCn))*lengthP;  xCsn2 = ones(size(zCn))*lengthP *(- cosd(60)); xCsn3 = ones(size(zCn))*lengthP *(- cosd(60));
yCsn1 = zeros(size(zCn));  yCsn2 = ones(size(zCn))*lengthP *(sind(60)); yCsn3 = ones(size(zCn))*lengthP *(- sind(60));
nodesCs1 = [xCsn1 yCsn1 zCsn1]; nodesCs2 =[xCsn2 yCsn2 zCsn2]; nodesCs3 =[xCsn3 yCsn3 zCsn3];
nodesCs = [nodesCs1;nodesCs2;nodesCs3];
NnodesCs1 = size(nodesC,1); NnodesCs2 = size(nodesC,1); NnodesCs3 = size(nodesC,1); 
% pontoon1
xP1n = xP1 - 0.5*(lengthP/NsectionP)*ones(size(xP1));   
xP1n = xP1n(2:end);
yP1n = zeros(size(xP1n));                
zP1n = Zs*ones(size(xP1n));
nodesP1 = [xP1n yP1n zP1n];   NnodesP1 = size(nodesP1,1);
% pontoon2 and 3
xP2n = -cosd(60)*xP1n;                    
yP2n = cosd(30)*xP1n;                     
zP2n = Zs*ones(size(yP2n));               
xP3n = xP2n;                              
yP3n = -yP2n;                             
zP3n = zP2n;
nodesP2 = [xP2n yP2n zP2n];   NnodesP2 = size(nodesP2,1);
nodesP3 = [xP3n yP3n zP3n];   NnodesP3 = size(nodesP3,1);
nodes = [nodesC; nodesCs1; nodesCs2; nodesCs3; nodesP1; nodesP2; nodesP3];  
Nnodes = [length(nodes), NnodesC, NnodesCs1, NnodesCs2, NnodesCs3, NnodesP1, NnodesP2, NnodesP3] ;     

%---------objectives-------------------
% nodesC, nodesCs1, nodesCs2, nodesCs3
% nodesP1, nodesP2, nodesP3
% nodes
% ------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% section 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate SIMO body input for A, B and F from radiation and diffraction pressures
% Now we defined 7 substructures, 4 columns and 3 pontoons
Nbody =7;
AB = cell(Nbody,1);
F = cell(Nbody,1);
Buoy = cell(Nbody,1);
for i =1:Nbody
    AB{i} = zeros(Npers,Nnodes(i+1),6,6);
    F{i} = zeros(Npers,Nbetas,Nnodes(i+1),6);
    Buoy{i} = zeros(Npers,Nnodes(i+1));
end
Lm = [L^2 L^2 L^2 L^3 L^3 L^3];                 % m for force. n and k for pressure and radiation coefficient.
Lkn = [L^-3 L^-3 L^-3 L^-4 L^-4 L^-4];          % for radiation. L k-n, see 6.18 equation. correct.     
count = 0;
n = zeros(Nbody,1);
for iper = 1:Npers   % number of periods , in total 60 
    n = zeros(Nbody,1);
    disp(['Run calculation period ' num2str(pers(iper))])  
    % loop, panels and separate. 2 seciton node. 3 beta for force
    for ipnl = 1:Npanels                    % loop through panels, all panels. this is panel number
        count = count + 1;                  % extract 5p file data. see below
        norm = [nx(ipnl), ny(ipnl), nz(ipnl)]; % normal direction
        area = AREA(ipnl);                     % panel area
       % Classify panels into different substructures and sections 
       %----------------------------------------------------------------------
        if any(pnl_C==ipnl)        % central column sequence
            [F{1}, AB{1}, Buoy{1}, n(1)] = Integration_Column(F{1},AB{1},Buoy{1},Nnodes(2),ipnl,ZCT,zC,dat5p,norm,n(1),count,Lkn,area,g,omega,iper,Nbetas,Npanels,Lm);
        elseif any(pnl_Cs1==ipnl)  % side column 1 sequence 
            [F{2}, AB{2}, Buoy{2}, n(2)] = Integration_Column(F{2},AB{2},Buoy{2},Nnodes(3),ipnl,ZCT,zC,dat5p,norm,n(2),count,Lkn,area,g,omega,iper,Nbetas,Npanels,Lm);
        elseif any(pnl_Cs2==ipnl)  % side column 2 sequence 
            [F{3}, AB{3}, Buoy{3}, n(3)] = Integration_Column(F{3},AB{3},Buoy{3},Nnodes(4),ipnl,ZCT,zC,dat5p,norm,n(3),count,Lkn,area,g,omega,iper,Nbetas,Npanels,Lm);
        elseif any(pnl_Cs3==ipnl)  % side column 3 sequence 
            [F{4}, AB{4}, Buoy{4}, n(4)] = Integration_Column(F{4},AB{4},Buoy{4},Nnodes(5),ipnl,ZCT,zC,dat5p,norm,n(4),count,Lkn,area,g,omega,iper,Nbetas,Npanels,Lm);
        elseif any(pnl_P1==ipnl)
            [F{5}, AB{5}, Buoy{5}, n(5)] = Integration_Pontoon(F{5},AB{5},Buoy{5},Nnodes(6),ipnl,XCT,YCT,ZCT,xP1,yP1,dat5p,norm,n(5),count,Lkn,area,g,omega,iper,Nbetas,Npanels,Lm);
        elseif any(pnl_P2==ipnl) 
            [F{6}, AB{6}, Buoy{6}, n(6)] = Integration_Pontoon(F{6},AB{6},Buoy{6},Nnodes(7),ipnl,XCT,YCT,ZCT,xP2,yP2,dat5p,norm,n(6),count,Lkn,area,g,omega,iper,Nbetas,Npanels,Lm);
        elseif any(pnl_P3==ipnl) 
            [F{7}, AB{7}, Buoy{7}, n(7)] = Integration_Pontoon(F{7},AB{7},Buoy{7},Nnodes(8),ipnl,XCT,YCT,ZCT,xP3,yP3,dat5p,norm,n(7),count,Lkn,area,g,omega,iper,Nbetas,Npanels,Lm); 
        end
       %-----------------------------------------------------------------------
      
    end 
     % each pers, have these panels
    count = count + Npanels*Nbetas; 
end

%% Set small numbers to 0
for i = 1:7
    F{i}(abs(real(F{i}))<1e-5 & abs(imag(F{i}))<1e-5) = 0;
    AB{i}(abs(real(AB{i}))<1e-5 & abs(imag(AB{i}))<1e-5) = 0;
end
   

%% section 4-Write results to files
% these result are wamit format, non-dimensional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% body 1 --> central column
% body 2,3,4 --> side columns
% body 5,6,7 --> pontoons
bodyname = ["Cc",'Cs1','Cs2','Cs3','P1','P2','P3'];
judge = input('if write files for SIMA ? -->');
if judge == 1
    Ncount = -Nnodes(2);
    for ibo = 1:7
        Ncount = Ncount + Nnodes(i+1);
        writefile(char(bodyname(ibo)), nodes(i), Nnodes(i+1), Ncount, folderout, pers, betas, Npers, Nbetas, mode_out, F{i}, AB{i})
    end
end
clear judge
%% section 5-check if the integration is correct

% Read the files with totals
% -------------------------------------------------------------------------
ABfile = [folder prefix '.1'];        % total AB come from the Wamit
dat = load(ABfile);                   
A_tot_WAMIT = zeros(Npers,6,6);
B_tot_WAMIT = zeros(Npers,6,6);
A = 1;
for kk = 1:length(dat(:,1))
   per = dat(kk,1);
   w = 2*pi/per;
   ii = dat(kk,2);
   jj = dat(kk,3);
   a = dat(kk,4);
   b = dat(kk,5);
   k = 3;
   if (ii>3); k = k+1; end
   if jj>3; k = k+1; end
   
   iper = find(pers==per);
   A_tot_WAMIT(iper,ii,jj) = a*rho*L^k;
   B_tot_WAMIT(iper,ii,jj) = b*w*rho*L^k;
end

Ffile= [folder prefix '.2'];
dat = load(Ffile);
F_tot_WAMIT = zeros(Npers,Nbetas,6);

for kk = 1:length(dat)
   per = dat(kk,1);
   beta = dat(kk,2);
   ii = dat(kk,3);
   m = 2;
   if ii>3; m=3;end
   f = (dat(kk,6)+dat(kk,7)*1i)*rho*g*A*L^m; 
   iper = find(pers==per);
   ibet = find(betas==beta);
   F_tot_WAMIT(iper,ibet,ii) = f; 
end
% -------------------------------------------------------------------------
%% 
% Sum up all the sections of all the bodies
AB_tot = zeros(Npers,6,6);
F_tot = zeros(Npers,Nbetas,6);
Buoy_tot = zeros(Npers,1);


for ibod =1:length(Nnodes)-1
    AB_temp = AB{ibod};
    F_temp = F{ibod};
    Buoy_temp = Buoy{ibod};
    [nodex, nodey, nodez] =  node_body (nodesC, nodesCs1, nodesCs2, nodesCs3, nodesP1, nodesP2, nodesP3, ibod);

    %loop over periods for each body

    for iper = 1:Npers
        for inod=1:Nnodes(ibod+1)
            
              % Buoyancy
              Buoy_tot(iper) = Buoy_tot(iper) + Buoy_temp(iper,inod);
             
              % translation terms
              for mode1=1:3 
                  for mode2 = 1:3
                    AB_tot(iper,mode1,mode2) = AB_tot(iper,mode1,mode2) + AB_temp(iper,inod,mode1,mode2)*rho*L^3; 
                  end
                  for ibeta = 1:Nbetas
                    F_tot(iper,ibeta,mode1) = F_tot(iper,ibeta,mode1) + F_temp(iper,ibeta,inod,mode1)*rho*g*A*L^2; 
                  end
              end
              
              for ibeta=1: Nbetas
              % forces in rotation
                F_tot(iper,ibeta,4) = F_tot(iper,ibeta,4)- F_temp(iper,ibeta,inod,2)*rho*g*A*L^m*nodez(inod) + F_temp(iper,ibeta,inod,3)*rho*g*A*L^m*nodey(inod); 
                F_tot(iper,ibeta,5) = F_tot(iper,ibeta,5)+ F_temp(iper,ibeta,inod,1)*rho*g*A*L^m*nodez(inod) - F_temp(iper,ibeta,inod,3)*rho*g*A*L^m*nodex(inod); 
                F_tot(iper,ibeta,6) = F_tot(iper,ibeta,6)+ F_temp(iper,ibeta,inod,2)*rho*g*A*L^m*nodex(inod) - F_temp(iper,ibeta,inod,1)*rho*g*A*L^m*nodey(inod); 
              end
              % rotation/translation
              AB_tot(iper,5,1) = AB_tot(iper,5,1) + AB_temp(iper,inod,1,1)*rho*L^3*nodez(inod) - AB_temp(iper,inod,1,3)*rho*L^3*nodex(inod); 
              AB_tot(iper,4,2) = AB_tot(iper,4,2) - AB_temp(iper,inod,2,2)*rho*L^3*nodez(inod) + AB_temp(iper,inod,2,3)*rho*L^3*nodey(inod); 
              AB_tot(iper,6,2) = AB_tot(iper,6,2) + AB_temp(iper,inod,2,2)*rho*L^3*nodex(inod) - AB_temp(iper,inod,1,2)*rho*L^3*nodey(inod); 
              % rotation
              AB_tot(iper,5,5) = AB_tot(iper,5,5) + AB_temp(iper,inod,1,1)*rho*L^3*(nodez(inod))^2 - 2*AB_temp(iper,inod,1,3)*rho*L^3*nodex(inod)*nodez(inod) + AB_temp(iper,inod,3,3)*rho*L^3*(nodex(inod))^2;
              AB_tot(iper,4,4) = AB_tot(iper,4,4) + AB_temp(iper,inod,2,2)*rho*L^3*(nodez(inod))^2 - 2*AB_temp(iper,inod,3,2)*rho*L^3*nodey(inod)*nodez(inod) + AB_temp(iper,inod,3,3)*rho*L^3*(nodey(inod))^2;
              AB_tot(iper,6,6) = AB_tot(iper,6,6) + AB_temp(iper,inod,1,1)*rho*L^3*(nodey(inod))^2 - 2*AB_temp(iper,inod,1,2)*rho*L^3*nodex(inod)*nodey(inod)+ AB_temp(iper,inod,2,2)*rho*L^3*(nodex(inod))^2;

        end
    end
    clear AB_temp F_temp
end

AB_tot(:,1,5) = AB_tot(:,5,1);
AB_tot(:,2,4) = AB_tot(:,4,2);
AB_tot(:,2,6) = AB_tot(:,6,2);

%% Plot checking
% excitation
plot_excitation(pers,F_tot_WAMIT, F_tot);


%% 
% selected added mass and damping plots (total)
plot_radiation(pers, A_tot_WAMIT,B_tot_WAMIT, AB_tot);

%% plot cross sections, RIFLEX Beam elements, SIMO bodies
close all
crd_Cs = [crd_Cs1;crd_Cs2;crd_Cs3];
crd_Pon = [crd_P1;crd_P2;crd_P3];
plot_model(idl,crd_C, crd_Cs, zC, lengthP, crd_Pon,crd_P1, crd_P2, crd_P3, xP1, XCT, YCT, ZCT, nodes,...
            Secs,SecP1, SecP2, SecP3, SecC, SecCs1, SecCs2, SecCs3);
%% write total result of multi-body model
%             dl = 1	 dl = 2	   dl = 3	 dl = 4   dl = 5
% Column x4	    16	        12	      8	       4         2
% Pontoon x3	16	        12	      8	       4         2
% Total	        108	        84	      56	   28        14

% mkdir multibody_total
folder_multi = 'multibody_total\';
write_multi_total(folder_multi, idl, Npers, pers, Nbetas, betas, F_tot, AB_tot)

%% Read from the file written with the integrated added mass of different multibodies
ifinish = input('if finished ? -->');

Npers = 60;
if ifinish == 1
    AB_tot_dl = zeros(Ndl, Npers, 6, 6);
    F_tot_dl = zeros(Ndl, Npers, Nbetas, 6);
     for idl = 1:Ndl
         AB_tot_file = [folder_multi 'AB_' num2str(idl) '.1'];
         F_tot_file = [folder_multi 'F_' num2str(idl) '.3'];1
         dat = load(AB_tot_file);
         dat2 = load(F_tot_file);
         for iper = 1:Npers
         AB_tot_dl(idl,iper,:,:) = reshape(dat(iper*36-35:iper*36,4)+1i*dat(iper*36-35:iper*36,5),[6,6])';
         F_tot_dl(idl,iper,:,:) = reshape(dat2((iper-1)*Nbetas*6+1:iper*Nbetas*6, 6) +...
                                          1i*dat2((iper-1)*Nbetas*6+1:iper*Nbetas*6, 7),[Nbetas,6]);           
         end
     end
     
     plot_div_mode(pers, A_tot_WAMIT,B_tot_WAMIT, AB_tot_dl, F_tot_WAMIT, F_tot_dl, Ndl, folder_multi)
end
%% Gravity and Buoyancy
clc
%thickness [m]
thk083 = 0.083;
thk12 = 0.12;
thk17 = 0.17;
thk20 = 0.20;
thk25 = 0.25;

% Mass ------------taken from IEA report
% unit[kg]
M_blades = 6.525e4*3;    % V 
M_hub = 1.9e5;         % V
M_naccelle = 6.3089e5; % V
M_RNA = 9.91e5;
M_tower = 1.263e6;     % V
M_plat= 1.7839e7;

M_hub + M_naccelle + M_blades;

M_tot = M_RNA + M_tower + M_plat;
M_report = 20093000;
if M_tot==M_report
    disp('total mass correct')
end
G_tot = M_tot * g;


%-------------------------------
M_steel = 3.914e6;
MBall_ironcon = 2.54e6;
MBall_sw = 1.13e7;
MBall = MBall_ironcon + MBall_sw;
M_plat_check = M_steel + MBall ;

% 33206*28 + 20288*28*3 + 23388*51.75*3
rho_iron = 7874;
Mcols = (pi*radC_s*2 *15*thk12 + pi*radC_s*2 *14*thk17 + pi*radC_s^2*thk25) * rho_iron;
Mcolc = pi*radC_c*2 *29  * thk083 * rho_iron;
Mpon = (M_plat - Mcols*3 - Mcolc)/3;

disp('Mass:')
disp('pontoon:-----------  central:------------   side:-----------')
disp([Mpon, Mcolc, Mcols])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Mass of Platform:')
disp([3*Mpon+3*Mcols+Mcolc])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
M_stiff = 3 * (0.91/2)^2*pi*lengthP*7085;
Mcoef_pon =  Mpon/(radC_s+lengthP);
Mcoef_colc = Mcolc/29;
Mcoef_cols = Mcols/28;
disp('Mcoeff:')
disp('pontoon:-----------  central:------------   side:-----------')
disp([Mcoef_pon,Mcoef_colc,Mcoef_cols])
%---------------------------------------











