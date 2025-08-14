function [F, AB, Buoy, n] = Integration_Column(Fin, ABin, Buoyin, Nnodes, ipnl, ZCT, zC, ...
                               dat5p, norm, n, count, Lkn, AREA, g, ...
                               omega, iper, Nbetas, Npanels, Lm)
% 
%
%
%
%
 F = Fin;
 AB = ABin;
 Buoy = Buoyin;
 heightC = 13;
 rho = 1025;
    for inode = 1:Nnodes          
    % from low to up, both the upper column and lower column
        
       if ZCT(ipnl) > zC(inode) && ZCT(ipnl) <= zC(inode+1)     % panel center between two node of section
          %Calculate added mass and damping matrix from radiation presssure 
          %p = pressure(dat5p,count);    
          p1 = dat5p(count,4)+1i*dat5p(count,5);                
          p2 = dat5p(count,6)+1i*dat5p(count,7);
          p3 = dat5p(count,8)+1i*dat5p(count,9);             
          p = [p1 p2 p3 0 0 0]; 
          clear p1 p2 p3
          rn = [0, 0, 0];                 % Set the rotational forces to 0.
          nrn = ([norm rn].*Lkn)'; 
          
          % 1. Added mass and Damping
          ABpan = nrn*p*AREA*g/(omega(iper))^2;          
          AB(iper,inode,:,:) = squeeze(AB(iper,inode,:,:)) + ABpan;  
          n = n+1;                     
          
          % 2. Exicatation of different wave directions
          for ibet = 1:Nbetas                                   
              f = dat5p(count+Npanels*ibet,5) + 1i*dat5p(count+Npanels*ibet,6);  
              Fpan = AREA./Lm.*[norm rn]*f; 
              clear f
              F(iper,ibet,inode,1) = Fpan(1) + F(iper,ibet,inode,1);  %Fx  
              F(iper,ibet,inode,2) = Fpan(2) + F(iper,ibet,inode,2);  %Fy
              F(iper,ibet,inode,3) = Fpan(3) + F(iper,ibet,inode,3);  %Fz
              F(iper,ibet,inode,4) = 0;
              F(iper,ibet,inode,5) = 0;
              F(iper,ibet,inode,6) = 0;
          end
              
           % 3. Buoy contributed by vertical pressure
           Buoy(iper,inode) =  Buoy(iper,inode) + norm(3) * AREA * rho * g * (-ZCT(ipnl));
       
        
%        %------------------------------------------------------------------------------------------------- 
%        elseif ZCT(ipnl) < -heightC % include the BN surface into the lowest section of the central column
%        %Calculate added mass and damping matrix from radiation presssure 
%           %p = pressure(dat5p,count);    
%           p1 = dat5p(count,4)+1i*dat5p(count,5);                
%           p2 = dat5p(count,6)+1i*dat5p(count,7);
%           p3 = dat5p(count,8)+1i*dat5p(count,9);             
%           p = [p1 p2 p3 0 0 0]; 
%           clear p1 p2 p3
%           rn = [0, 0, 0];                 % Set the rotational forces to 0.
%           nrn = ([norm rn].*Lkn)'; 
%           
%           % 1. Added mass and Damping
%           ABpan = nrn*p*AREA*g/(omega(iper))^2;          
%           AB(iper,1,:,:) = squeeze(AB(iper,1,:,:)) + ABpan;  
%           n = n+1;                       
%           
%           % 2. Exicatation of different wave directions
%           for ibet = 1:Nbetas                                   
%               f = dat5p(count+Npanels*ibet,5) + 1i*dat5p(count+Npanels*ibet,6);  
%               Fpan = AREA./Lm.*[norm rn]*f; 
%               clear f
%               F(iper,ibet,1,1) = Fpan(1) + F(iper,ibet,1,1);  %Fx  
%               F(iper,ibet,1,2) = Fpan(2) + F(iper,ibet,1,2);  %Fy
%               F(iper,ibet,1,3) = Fpan(3) + F(iper,ibet,1,3);  %Fz
%               F(iper,ibet,1,4) = 0;
%               F(iper,ibet,1,5) = 0;
%               F(iper,ibet,1,6) = 0;
%           end
%           % 3. Buoy contributed by vertical pressure
%           Buoy(iper,1) =  Buoy(iper,1) + norm(3) * AREA * rho * g * (-ZCT(ipnl));
          
       end
      
    end
    
    
