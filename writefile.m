function writefile(bodyname, nodes, Nnodes, Ncount, folderout, pers, betas, Npers, Nbetas, mode_out,...
                    F, AB)

for body = 1:Nnodes               
    %Nbody = body + 0;                % name, total number 
    disp([bodyname ' ' num2str(body) ]) % only for display
   
    countdiff = 0;
    countrad = 0;
    fileID = fopen([folderout 'body' num2str(body+Ncount) '.1'],'w');  
    fileID2 = fopen([folderout 'body' num2str(body+Ncount) '.3'],'w');
    disp(['In body.retf:' '[' num2str(body+Ncount) ']']);
    for iper2 = 1:Npers
        for ibeta = 1:Nbetas
            for mode = 1:mode_out

                % Print excitation forces to file
                countdiff = countdiff + 1;
                Mod = abs(F(iper2,ibeta,body,mode));            
                Pha = angle(F(iper2,ibeta,body,mode))*180/pi;
                Re = real(F(iper2,ibeta,body,mode));
                Im = imag(F(iper2,ibeta,body,mode));
                fprintf(fileID2,'%.4e\t%.4e\t%i\t%.4e\t%.4e\t%.4e\t%.4e\n',[pers(iper2), betas(ibeta), mode, Mod, Pha, Re, Im]);
                    
                % Print added mass and damping to file
                countrad = countrad + 1;
                if ibeta==1                
                   for j = 1:mode_out
                       %if abs(AB_C(iper2,body,mode,j)) > eps      
                          Aij = real(AB(iper2,body,mode,j));     
                          Bij = -imag(AB(iper2,body,mode,j));      
                          fprintf(fileID,'%.4e\t%i\t%i\t%.4e\t%.4e\n',[ pers(iper2), mode, j, Aij, Bij ]);
                       %end
                   end
                end
                clear Mod Pha Re Im Aij Bij
            end
        end
    end

    fclose(fileID);
    fclose(fileID2);

%     copyfile([folderout 'semi.mmx'],[folderout 'body' num2str(body) '.mmx'])    
%     
%     %Replace all spar with bodyXX in file to get more redable simo body names
%     fidi  = fopen([folderout 'spar.out'],'r');        
%     f=fread(fidi,'*char')';                            
%     fclose(fidi);  
%     f = strrep(f,'semi',['body' num2str(body)]);      
%     fidi  = fopen([folderout 'Cbody' num2str(body) '.out'],'w');   
%     fprintf(fidi,'%s',f);       
%     fclose(fidi);

end