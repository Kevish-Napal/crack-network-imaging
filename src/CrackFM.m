function ICrack = CrackFM(ListX,ListY,iter_k,Mreg)
% Computes the indicator of the crack defined in the Paper Haddar-Yosra

global Feta HEST List_k


[nbcapteur,~,~] = size(Feta);
Feta = 2*pi/nbcapteur*Feta;
capteur = 0:2*pi/nbcapteur:2*pi-2*pi/nbcapteur;

k = List_k(iter_k);

[ZX,ZY] = meshgrid(ListX,ListY);


% Calcul de Fsharp
[Vr,Dr]=eig((Feta(:,:,iter_k)+Feta(:,:,iter_k)')/2);
[Vi,Di]=eig((Feta(:,:,iter_k)-Feta(:,:,iter_k)')/(2*sqrt(-1)));
Fsharp=Vr*abs(Dr)*Vr^-1+Vi*abs(Di)*Vi^-1;
[V,D] = eig(Fsharp);
lambda = abs(diag(D));% these are positive reals, abs is just is here to erase the small im part

ICrack = zeros(size(ZX));

dx=cos(capteur);dy=sin(capteur);               % coordonnées cartésiennes des capteurs.
NbPoints = numel(ZX);

zeta = -1:.001:1;
zeta2 = zeta.^2;
        
    parfor i = 1:NbPoints 
   
     %disp([num2str(i) '/' num2str(NbPoints)]);
        z = [ZX(i) ZY(i)];
        
        eta_d = 1/sqrt(8.*pi*k)*exp(1i*pi/4.); 
        L = 0.001;
        
        % computation of Phi_{z,L}^infty defined in yosra just before
        % equation (37). rhs corresponds to alpha = 0, rhs1 to beta = 0 and
        % nu = e1, rhs2 to beta = 0 and nu = e2.
        rhs=transpose(eta_d*L*(exp(-sqrt(-1).*k.*(dx.*z(1)+dy.*z(2))))); 
        rhs1 = -1i * k * dx' .* rhs;
        rhs2 = -1i * k * dy' .* rhs;
        
        
        % First Picard Criteria (norm of g)
        sc=V'*rhs;
        GUP=conj(sc).*sc;
        Gl = GUP./lambda;
        Ng = sum(Gl(1:Mreg));
        
        % gnu1 and gnu2
        sc=V'*rhs1;
        gnu1 = V*(sc./sqrt(lambda));
        gnu1 = gnu1(1:Mreg);
        Ngnu1 = dot(gnu1,gnu1);
        
        % gnu1 and gnu2
        sc=V'*rhs2;
        gnu2 = V*(sc./sqrt(lambda));
        gnu2 = gnu2(1:Mreg);
        Ngnu2 = dot(gnu2,gnu2);
        
        gnu1dotgnu2 = dot(gnu1,gnu2);
        
        Ngnu = min( Ngnu1*zeta2 + Ngnu2*(1-zeta2) + 2*real(gnu1dotgnu2)*zeta.*sqrt(1-zeta2) );

        
        ICrack(i) = 1/sqrt(Ng) + 1/sqrt(Ngnu);
        
     
        
    end
    
    

end
