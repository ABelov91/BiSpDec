function [z,h,pair_num,Nums] = geometry
% d1 - толщина первого слоя в периоде, d2 - второго и т.д.
% d_substr - толщина слоя-подложки
% d_cover - толщина слоя-покрытия
% d_front, d_back - слои до и после структуры
% pair_num - число пар слоев

global problem;

global ppSIO2;
global ppTA2O5;

global theta_inc;
global omega_mid;

global thick1;
global thick2;

global fluct;

if( problem == 0 )
    n_front = 1;
    n_cover = 1;
    n1 = 2;
    n2 = 4;
    n_substr = 1;
elseif( problem == 1 )
    n_front = 1;
    n_cover = 1;
    n1 = real( ppval(ppSIO2,omega_mid)^0.5 );
    n2 = real( ppval(ppTA2O5,omega_mid)^0.5 );
    n_substr = 1;
elseif( problem == 2 || problem == 3 || problem == 4 )
    n_front = 1.49;
    n_cover = 1.49;
    n1 = real( ppval(ppSIO2,omega_mid)^0.5 );
    n2 = real( ppval(ppTA2O5,omega_mid)^0.5 );
    n_substr = real( ppval(ppSIO2,omega_mid)^0.5 );
elseif( problem == 5 )
    n_front = 1;
    n_cover = real( ppval(ppSIO2,omega_mid)^0.5 );
    n1 = real( ppval(ppTA2O5,omega_mid)^0.5 );
    n2 = real( ppval(ppSIO2,omega_mid)^0.5 );
    n_substr = 1;
end

if( problem == 0 )
    d_front = 0.05;
    d_cover = 0.05;
    pair_num = 7;
    d1 = 1/8;
    d2 = 1/16;
    d_substr = 0.1;
    d_back = 0.05;
    
    N_front = 1;
    N_cover = 1;
    N1 = 1;
    N2 = 1;
    N_substr = 1;
    N_back = 1;
elseif( problem == 1 )
    d_front = 0.025;
    d_cover = 0.025;
    pair_num = 7;
    d1 = 0.130;
    d2 = 0.092;
    d_substr = 0.260;
    d_back = 0.1;
    
    N_front = 1;
    N_cover = 5;
    N1 = 2;
    N2 = 2;
    N_substr = 1;
    N_back = 1;
elseif( problem == 2 )
    d_front = 0.1;
    d_cover = 0.1;
    pair_num = 7;
    d1 = 0.160;
    d2 = 0.112;
    d_substr = 0.260;
    d_back = 0.8;
    
    N_front = 1;
    N_cover = 1;
    N1 = ceil(d1*n1/d_front);
    N2 = ceil(d2*n2/d_front);
    N_substr = ceil(d_substr*n_substr/d_front);
    N_back = 8;

elseif( problem == 3 )
    % толщины слоев
    d_front = 0.1;
    d_cover = 0.1;
    pair_num = 7;
    d1 = 0.207;
    d2 = 0.145;
    d_substr = 0.185;
    d_back = 0.8;
    
    % числа шагов в каждом слое
    N_front = 1;
    N_cover = 1;
    N1 = ceil(d1*n1/d_front);
    N2 = ceil(d2*n2/d_front);
    N_substr = ceil(d_substr*n_substr/d_front);
    N_back = 8;
    
elseif( problem == 4 )
    d_front = 0.1;
    d_cover = 0.1;
    pair_num = 7;
    d1 = thick1;
    d2 = thick2;
    d_substr = 0.260;
    d_back = 0.8;
    
    N_front = 1;
    N_cover = 1;
    N1 = ceil(d1*n1/d_front);
    N2 = ceil(d2*n2/d_front);
    N_substr = ceil(d_substr*n_substr/d_front);
    N_back = 8;

elseif( problem == 5 )
    d_front = 0.1;
    d_cover = 0.260;%thick2;
    pair_num = 7;
    d1 = thick1;
    d2 = thick2;
    d_substr = 0.1;
    d_back = 0.1;
    
    N_front = 1;
    N_cover = ceil(d_cover*n_cover/d_front);
    N1 = ceil(d1*n1/d_front);
    N2 = ceil(d2*n2/d_front);
    N_substr = 1;
    N_back = 1;
end

sin_theta_cover = n_front/n_cover*sin(theta_inc);
cos_theta_cover = sqrt( 1 - sin_theta_cover^2 );

sin_theta1 = n_front/n1*sin(theta_inc);
cos_theta1 = sqrt( 1 - sin_theta1^2 );

sin_theta2 = n_front/n2*sin(theta_inc);
cos_theta2 = sqrt( 1 - sin_theta2^2 );

sin_theta_substr = n_front/n_substr*sin(theta_inc);
cos_theta_substr = sqrt( 1 - sin_theta_substr^2 );

d_cover = d_cover*abs(cos_theta_cover);
d1 = d1*abs(cos_theta1);
d2 = d2*abs(cos_theta2);
d_substr = d_substr*abs(cos_theta_substr);

N0 = N_front + N_cover + pair_num*N1 + pair_num*N2 + N_substr + N_back;

h = zeros(1,N0);
for n = 1:N_front
    h(n) = d_front/N_front;
end
for n = 1:N_cover
    h(N_front + n) = d_cover/N_cover;
end
for m = 1:pair_num
    count = N_front + N_cover + (m-1)*(N1+N2);
    for n = 1:N1
        h(count + n) = (d1 + fluct(2*m-1))/N1;
    end
    count = N_front + N_cover + (m-1)*(N1+N2) + N1;
    for n = 1:N2
        h(count + n) = (d2 + fluct(2*m))/N2;
    end
end
count = N_front + N_cover + pair_num*(N1+N2);
for n = 1:N_substr
    h(count + n) = (d_substr + fluct(2*pair_num+1))/N_substr;
end
count = N_front + N_cover + pair_num*(N1+N2) + N_substr;
for n = 1:N_back
    h(count + n) = d_back/N_back;
end

z = zeros(1,N0+1);
for n = 1:N0
    z(n+1) = z(n) + h(n);
end

Nums = [N_front, N_cover, N1, N2, N_substr, N_back];

end