function eps = diel(Z, omega)
% eps1 - eps первого слоя в периоде, eps2 - второго и т.д.
% eps_substr - eps слоя-подложки
% eps_front, eps_back - слои до и после структуры
% pair_num - число пар слоев

global problem;
global ppSIO2;
global ppTA2O5;

if( problem == 0 )
    eps_front = 1;
    eps_cover = 1;
    eps1 = 4;
    eps2 = 16;
    eps_substr = 1;
    eps_back = 1;
elseif( problem == 1 )
    eps_front = 1;
    eps_cover = 1;
    eps1 = ppval(ppSIO2,omega);
    eps2 = ppval(ppTA2O5,omega);
    eps_substr = 1;
    eps_back = 1;
elseif( problem == 2 || problem == 3 || problem == 4 )
    eps_front = 1.49^2;
    eps_cover = 1.49^2;
    eps1 = ppval(ppSIO2,omega);
    eps2 = ppval(ppTA2O5,omega);
    eps_substr = ppval(ppSIO2,omega);
    eps_back = 1;
elseif( problem == 5 )
    eps_front = 1;
    eps_cover = ppval(ppSIO2,omega);
    eps1 = ppval(ppTA2O5,omega);
    eps2 = ppval(ppSIO2,omega);
    eps_substr = 1;
    eps_back = 1;
end

[z,~,pair_num,N] = geometry;
N_front = N(1);
N_cover = N(2);
N1 = N(3); N2 = N(4);

if( Z <= z(1+N_front) )
    eps = eps_front;
end
if( z(1+N_front) <= Z && Z <= z(1+N_front+N_cover) )
    eps = eps_cover;
end
for m = 1:pair_num
    count = 1 + N_front + N_cover + (m-1)*(N1+N2);
    if( z(count) <= Z && Z <= z(count + N1) )
        eps = eps1;
    end
    if( z(count + N1) <= Z && Z <= z(count + N1 + N2) )
        eps = eps2;
    end
end
N_substr = N(5);
count = 1 + N_front + N_cover + pair_num*(N1+N2);
if( z( count ) <= Z && Z <= z( count + N_substr ) )
    eps = eps_substr;
end
if( Z >= z( count + N_substr ) )
    eps = eps_back;
end 

end