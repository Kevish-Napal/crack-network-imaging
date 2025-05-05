function [DE,pmax] = DirichletEigenvalues(EigMax)
% ListOfDirichletEigenvalues = DirichletEigenvalues(EigMax).
% Returns the ordered list of all Dirichlet eigenvalues 
% (of the Laplace operator on the unit ball) contained in interval [0,EigMax].
% The precision is eps =  2.2204e-16.
% pmax is the maximum integer p for which the Bessel function J_p(x) has a zero
% in the interval [0,EigMax].



%% Initialisation: zeros of x -> J_0(x)

    p=0;
    f = @(x) besselj(p,x);
    dx=0.01;   
    DEp = fMultiZeros(f,0,EigMax,dx);
    DE = DEp;

%% Calcul des zeros de x -> J_p+1(x) a partir des zeros de x -> J_p(x)
% Les zeros de J_p+1(x) sont intercalés entre les zeros de J_p(x) 
    while isempty(DEp)==0
        p = p+1;
        f = @(x) besselj(p,x);
        DEpp = zeros(1,length(DEp)-1);
        for j = 1:(length(DEp)-1)
            DEpp(j) = fzero(f,[DEp(j),DEp(j+1)]);
        end
     
        % il faut rajouter les zeros de droite a la main
        DEpp = [DEpp,fMultiZeros(f,DEp(end),EigMax,0.01)];
        DEp = DEpp;
        DE = [DE,DEp];     
    end

    DE = sort(DE);
    pmax = p-1;
    if pmax ==-1
        disp(['No Dirichlet Eigenvalues in interval [0,' num2str(EigMax) ']']);
    end


