function poly = DiskSwipe(BX,BY,Rin,Ntheta)
% Rend une image lisible de la matrice M construite dans le main
    
    poly = [];
    theta = 0:2*pi/Ntheta:2*pi*(1-1/Ntheta);
    
    for b = 1:numel(BX)
        
        polytemp = polyshape(Rin*cos(theta)+BX(b),Rin*sin(theta)+BY(b));
        poly = [poly polytemp];
        
    end   

end