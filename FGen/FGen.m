function F = FGen(k,Rin,a,NbCapteur,BC)
% F = FGen(k,Rin,a,NbCapteur,BC) Creates a discretization of the far field operator associated to the
% exterior obstacle problem to the disk of radius Rin of center a = [ax ay]
% and of chosen BoundaryCondition = [BC1 BC2](linear combination of Dirichlet and Neumann).

BC1 = BC(1);
BC2 = BC(2);

    JPrime = @(p,x) -0.5*(besselj(p+1,x) - besselj(p-1,x));
    HPrime = @(p,x) -0.5*(besselh(p+1,x) - besselh(p-1,x));
    xp = @(p,k,Rin) - (BC1*besselj(p,k*Rin) + BC2*k*JPrime(p,k*Rin)) ...
        /(BC1*besselh(p,k*Rin) + BC2*k*HPrime(p,k*Rin));



    ThetaCapteur = 2*pi/NbCapteur*(0:NbCapteur-1);

    uinfty = 0;

    p = 0;
    add = xp(p,k,Rin)*cos(p*ThetaCapteur);


      while (norm(add,'inf') > 1e-12)
        uinfty = uinfty+add;
        p   = p + 1;
        add = 2*xp(p,k,Rin)*cos(p*ThetaCapteur);
      end


    uinfty = uinfty * sqrt(2/(pi*k)) * exp(-1.i*pi/4);




    % we use the symetry of the centered disk to compute the far field
    % operator F. The modification of the far field in case of uncentered
    % disk is made in the next section.
    F = transpose(toeplitz([uinfty(1) fliplr(uinfty(2:end))],uinfty));


    %% Translation if required:

    if norm(a) > 0
        [d,xhat] = meshgrid(ThetaCapteur);
        [Xd,Yd] = pol2cart(d,1);
        [Xxhat,Yxhat] = pol2cart(xhat,1);
        sx = Xd - Xxhat; sy = Yd - Yxhat;
        T = exp(1.i*k*(a(1)*sx + a(2)*sy));
        F = T.*F;
    end

end
