function [IMG,xmin,xmax,dx,ymin,ymax,dy] = Affichage(BX,BY,M,Rin)
% Rend une image lisible de la matrice M construite dans le main
    xmin = -4 ; xmax = 4; 
    ymin = -4 ; ymax = 4;
    dx = Rin/((xmax - xmin)*10);
    dy = Rin/((ymax - ymin)*10);
    
    [X,Y] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);
    IMG = zeros(size(X));
    for b = 1:numel(M)
        
        bool = (X-BX(b)).^2 + (Y-BY(b)).^2 < Rin^2;
        IMG(bool) = M(b);
        
    end   

end