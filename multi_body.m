function multi_body(NsectionC, NsectionP, d1)

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

