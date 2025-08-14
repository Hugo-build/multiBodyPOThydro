function write_multi_total(folder_multi, idl,  Npers, pers, Nbetas, betas, F_tot, AB_tot)

% This function is to save the integrated result of multibody model with
% different division of sections
fileID = fopen([folder_multi 'AB_' num2str(idl) '.1'],'w'); 
fileID2 = fopen([folder_multi 'F_' num2str(idl) '.3'],'w');
for iper2 = 1:Npers
    for imode =1:6
        for ibeta = 1:Nbetas
            Mod = abs(F_tot(iper2,ibeta,imode));            
            Pha = angle(F_tot(iper2,ibeta,imode))*180/pi;
            Re = real(F_tot(iper2,ibeta,imode));
            Im = imag(F_tot(iper2,ibeta,imode));
            fprintf(fileID2,'%.4e\t%.4e\t%i\t%.4e\t%.4e\t%.4e\t%.4e\n',[pers(iper2), betas(ibeta), imode, Mod, Pha, Re, Im]);  
        end
        for jmode = 1:6
            fprintf(fileID,'%.4e\t%i\t%i\t%.4e\t%.4e\n', [ pers(iper2), imode, jmode, real(AB_tot(iper2,imode,jmode)), imag(AB_tot(iper2,imode,jmode)) ]);
        end
    end
end
