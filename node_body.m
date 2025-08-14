function [nodex, nodey, nodez] =  node_body (nodesC, nodesCs1, nodesCs2, nodesCs3, nodesP1, nodesP2, nodesP3, ibod)

% This is to return the corresponding coordinates of the sections of the body 
if ibod == 1
        nodex = nodesC(:,1); nodey = nodesC(:,2) ;nodez = nodesC(:,3);
    elseif ibod == 2 
        nodex = nodesCs1(:,1); nodey = nodesCs1(:,2) ;nodez = nodesCs1(:,3);
    elseif ibod == 3 
        nodex = nodesCs2(:,1); nodey = nodesCs2(:,2) ;nodez = nodesCs2(:,3);
    elseif ibod == 4  
        nodex = nodesCs3(:,1); nodey = nodesCs3(:,2) ;nodez = nodesCs3(:,3);
    elseif ibod == 5
        nodex = nodesP1(:,1); nodey = nodesP1(:,2) ;nodez = nodesP1(:,3);
    elseif ibod == 6
        nodex = nodesP2(:,1); nodey = nodesP2(:,2) ;nodez = nodesP2(:,3);
    elseif ibod == 7
        nodex = nodesP3(:,1); nodey = nodesP3(:,2) ;nodez = nodesP3(:,3);
end