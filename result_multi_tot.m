close all
clear all
clc
tic
set(0,'defaultfigurecolor','w')
format long

prefix = 'semi';              
folder = 'C:\Users\Fish\OneDrive - NTNU\NTNU_postgraduate\Marine Structure\Master_Thesis\wamit_integrate\Wam_semi\'; 

%% Read WAMIT output
datGDF = dlmread([folder prefix '.GDF'],'',1,0);    % geometric data file
datpnl = dlmread([folder prefix '.pnl']);           % panel data file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
dat1 = dlmread([folder prefix '.1']);               % added mass and damping coefficients
dat2 = dlmread([folder prefix '.2']);               % exciting forces from Haskind relations
dat3 = dlmread([folder prefix '.3']);               % exciting forces from diffraction potential
dat5p = dlmread([folder prefix '.5p']);             % Hydrodynamic pressure on body surface

pers = flipud(unique(dat5p(:,1)));             
Npers = length(pers);                          % number of periods
rho = 1025;
L = 1;

%% Read the files with totals directly from WAMIT
% ------------------------------------------------------------------------

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
%% Read from the file written with the integrated added mass of different multibodies
AB_tot_file = 
AB_tot = zeros(4, Npers, 6, 6);
 for idl = 1:4
     for iper = 1:Npers
     AB_tot(idl,iper, ,) = 
     end
 end







