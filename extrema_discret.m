function L_A = extrema_discret(w_so, M, N, n_oct, n_spo)
%extrema_discret détermine les extrema discrets du DoG
%   ARGUMENTS:
%   	w_so: DoG digital (o=1,...,n_oct et s=0,...,n_spo+2)
%   OUTPUT:
%       L_A: liste des extremas du DoG

%% Initialisation
L_A = [];

%% Boucle de comparaison
for o=1:n_oct
    Mo = floor(M/2^(o-1));
    No = floor(N/2^(o-1));
    for s=1:n_spo
       for m=1:Mo-2
          for n=1:No-2
            largeur = [-1,1];
            test1p = prod(prod(w_so{s+1,o+1}(m+1,n+1)>w_so{s+1,o+1}(m+1+largeur,n+1+largeur)));
            test1p = test1p*prod(w_so{s+1,o+1}(m+1,n+1)>w_so{s+1,o+1}(m+1,n+1+largeur));
            test1p = test1p*prod(w_so{s+1,o+1}(m+1,n+1)>w_so{s+1,o+1}(m+1+largeur,n+1));
            
            largeur = -1:1;
            test2p = prod(prod(w_so{s+1,o+1}(m+1,n+1)>w_so{s+2,o+1}(m+1+largeur,n+1+largeur)));
            test3p = prod(prod(w_so{s+1,o+1}(m+1,n+1)>w_so{s,o+1}(m+1+largeur,n+1+largeur)));
            testp = test1p*test2p*test3p;
            
            largeur = [-1,1];
            test1m = prod(prod(w_so{s+1,o+1}(m+1,n+1)<w_so{s+1,o+1}(m+1+largeur,n+1+largeur)));
            test1m = test1m*prod(w_so{s+1,o+1}(m+1,n+1)<w_so{s+1,o+1}(m+1,n+1+largeur));
            test1m = test1m*prod(w_so{s+1,o+1}(m+1,n+1)<w_so{s+1,o+1}(m+1+largeur,n+1));
            
            largeur = -1:1;
            test2m = prod(prod(w_so{s+1,o+1}(m+1,n+1)<w_so{s+2,o+1}(m+1+largeur,n+1+largeur)));
            test3m = prod(prod(w_so{s+1,o+1}(m+1,n+1)<w_so{s,o+1}(m+1+largeur,n+1+largeur)));
            testm = test1m*test2m*test3m;
            
            
            if (testp==1) || (testm==1)
                L_A = [L_A;o,s,m,n];
            end
          end
       end
    end
end

end

