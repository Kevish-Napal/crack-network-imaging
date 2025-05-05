% Ce code permet d'extraire les données produites par FreeFem++ contenue
% dans le code Sol.edp. Le fichier DnuRef.mat contient DnuRef = d_nu(w-Phi_z) où w est
% la solution de l'équation de Helmholtz avec données de Dirichlet égale à
% Phi_z sur le disque centré en l'origine de rayon R discrétisé en Nh points. List_k contient
% les samplings en k et List_z les samplings en z.

close all
clear all
dossier = '' ;%'/home/labcmap/napal/Documents/numerique/Contour01.2019/CrackVertical/';



RealTemp = dlmread([dossier 'DnuReal.dat']);
ImageTemp = dlmread([dossier 'DnuImage.dat']);

R = RealTemp(1,1);
Nh = RealTemp(1,5);

shift = (Nh + 5)/5; %nombre de ligne par itération
Extract_Real = zeros(shift,5);
Extract_Image = zeros(shift,5);

Nk = 1;
NbZ = 1;
c = 0;

for iterk = 1:Nk
    for iterz = 1:NbZ
        
        Extract_Real = RealTemp(1+c*shift:(c+1)*shift,:);
        Extract_Image = ImageTemp(1+c*shift:(c+1)*shift,:);
        
        List_k(iterk) = Extract_Real(1,2);
        List_z(1,iterz) = Extract_Real(1,3);
        List_z(2,iterz) = Extract_Real(1,4);
        
        DnuRef(:,iterz,iterk) = reshape(transpose(Extract_Real(2:end,:)+1.i*Extract_Image(2:end,:)),[],1);
        
        
        
        c = c+1;
   end
end
save([dossier 'DnuRef.mat'],'R','Nh','List_k','List_z','DnuRef')
