function [I_LSM,I_GLSM] = CrackLSM(ListX,ListY,iter_k)
% Computes the indicator of the crack defined in the Paper Haddar-Yosra

global Feta HEST List_k

[nbcapteur,~,~] = size(Feta);
capteur = 0:2*pi/nbcapteur:2*pi-2*pi/nbcapteur;

h = HEST(iter_k);
k = List_k(iter_k);

[ZX,ZY] = meshgrid(ListX,ListY);


% Calcul de Fsharp^(1/2)
[Vr,Dr]=eig((Feta(:,:,iter_k)+Feta(:,:,iter_k)')/2);
[Vi,Di]=eig((Feta(:,:,iter_k)-Feta(:,:,iter_k)')/(2*sqrt(-1)));
F = Vr*sqrt(abs(Dr))*Vr^-1+Vi*sqrt(abs(Di))*Vi^-1;

[Ueta,Seta,Veta]=svd(F);             % Ueta et Veta sont des matrices de rotations, Seta une matrice diagonale contenant les valeurs singulières de Feta.
SIGMAN=diag(Seta);                      % SIGMAN est un vecteur contenant la diagonale de Seta.
SQ=SIGMAN.^2;                           % Vecteur contenant les valeurs propres de (Feta)'*Feta.
Us=Ueta';
options=optimset('Display','off');   

ICrack = zeros(size(ZX));

dx=cos(capteur);dy=sin(capteur);               % coordonnées cartésiennes des capteurs.
NbPoints = numel(ZX);

zeta = -1:.001:1;
zeta2 = zeta.^2;
        
    parfor i = 1:NbPoints 
   
     %disp([num2str(i) '/' num2str(NbPoints)]);
        z = [ZX(i) ZY(i)];
        
        eta_d = 1/sqrt(8.*pi*k)*exp(1i*pi/4.); 
        L = 1;
        
        % computation of Phi_{z,L}^infty defined in yosra just before
        % equation (37). rhs corresponds to alpha = 0, rhs1 to beta = 0 and
        % nu = e1, rhs2 to beta = 0 and nu = e2.
        rhs=transpose(eta_d*L*(exp(-sqrt(-1).*k.*(dx.*z(1)+dy.*z(2))))); 
        rhs1 = -1i * k * dx' .* rhs;
        rhs2 = -1i * k * dy' .* rhs;
        
        
        % First Picard Criteria (norm of g)
         % Definit la fonction dfun ici pour Morozov
         sc=Us*rhs;
        GUP=conj(sc).*sc;
        f = @(x) sum( GUP.*(x^2 - (h^2 * SQ))./((SQ+x).^2));
        
        
        
        gammamin=h*min(SIGMAN);    % borne inf pour la recherche de alpha* de morozov.
        gammamax=h*max(SIGMAN);    % le alpha(delta) de morozov se trouve dans l'intervalle [gammamin,gammamax].
      
        
        [a,fval,exitflag] = fzero(f,[gammamin,gammamax],options); %Faire par dichotomie pour voir (nombre de dichotomie?)
        % a =morozov,     fval = dfun(a),     exitflag < 0 => aucuns zero 

        if a<0||exitflag<0
            disp(['Problem: a= ',num2str(a),' Resetting a=0'])
            disp(['Problem: exitflag= ' num2str(exitflag)])
            a=0;
        end

      
        Valxi=(SIGMAN.*sc)./(SQ+a);
        g = Veta*Valxi;
        Ng = abs(dot(F*g,g))+h*dot(g,g);
        
        
        % gnu1 and gnu2
        sc=Us*rhs1;
         GUP=conj(sc).*sc;
        f = @(x) sum( GUP.*(x^2 - (h^2 * SQ))./((SQ+x).^2));
        
        
        
        gammamin=h*min(SIGMAN);    % borne inf pour la recherche de alpha* de morozov.
        gammamax=h*max(SIGMAN);    % le alpha(delta) de morozov se trouve dans l'intervalle [gammamin,gammamax].
      
        
        [a,fval,exitflag] = fzero(f,[gammamin,gammamax],options); %Faire par dichotomie pour voir (nombre de dichotomie?)
        % a =morozov,     fval = dfun(a),     exitflag < 0 => aucuns zero 

        if a<0||exitflag<0
            disp(['Problem: a= ',num2str(a),' Resetting a=0'])
            disp(['Problem: exitflag= ' num2str(exitflag)])
            a=0;
        end

      
        Valxi=(SIGMAN.*sc)./(SQ+a);
        gnu1 = Veta*Valxi;
        Ngnu1 = abs(dot(F*gnu1,gnu1))+h*dot(gnu1,gnu1);
        
        % gnu1 and gnu2
        sc=Us*rhs2;
         GUP=conj(sc).*sc;
        f = @(x) sum( GUP.*(x^2 - (h^2 * SQ))./((SQ+x).^2));
        
        
        
        gammamin=h*min(SIGMAN);    % borne inf pour la recherche de alpha* de morozov.
        gammamax=h*max(SIGMAN);    % le alpha(delta) de morozov se trouve dans l'intervalle [gammamin,gammamax].
      
        
        [a,fval,exitflag] = fzero(f,[gammamin,gammamax],options); %Faire par dichotomie pour voir (nombre de dichotomie?)
        % a =morozov,     fval = dfun(a),     exitflag < 0 => aucuns zero 

        if a<0||exitflag<0
            disp(['Problem: a= ',num2str(a),' Resetting a=0'])
            disp(['Problem: exitflag= ' num2str(exitflag)])
            a=0;
        end

      
        Valxi=(SIGMAN.*sc)./(SQ+a);
        gnu2 = Veta*Valxi;
        Ngnu2 = abs(dot(F*gnu2,gnu2))+h*dot(gnu2,gnu2);
        
        gnu1dotgnu2 = dot(F*gnu1,gnu2)+h*dot(gnu1,gnu2);
        
        Ngnu = min( Ngnu1*zeta2 + Ngnu2*(1-zeta2) + 2*real(gnu1dotgnu2)*zeta.*sqrt(1-zeta2) );

        
        I_LSM(i) = 1/sqrt(Ng) + 1/sqrt(Ngnu);
        
     
        
    end
    
    

end
