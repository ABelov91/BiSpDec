function [Fields, Fields_refl, reflection, reflection_norm, transmittance, cos_theta] = stationary_1d_single_mesh(param, z)
global theta_inc;
global polarization;

% param = [N1,N2, E0, I0, omega]
N1 = param(1);
N2 = param(2);
E0 = param(3);
I0 = param(4);
omega = param(5);

% параметры задачи

c = 2.997924580;
a = z(end);

N = length(z) - 1;

zmid = zeros(1,N);
h = zeros(1,N);
for n = 1:N
	h(n) = z(n+1) - z(n);
end
for n = 1:N
    zmid(n) = 0.5*z(n) + 0.5*z(n+1);
end

eps = zeros(1,N);
mu  = zeros(1,N);
j   = zeros(1,N);
j_surf = zeros(1,N+1);

for n = 1:N
    eps(n) = diel( zmid(n),omega );
    mu(n)  = magn( zmid(n),omega );
    j(n)   = curr( zmid(n) );    
end
for n = 1:N+1
    j_surf(n) = I0*surf_curr( z(n) );
end

cos_theta = cos(theta_inc)*ones(1,N);
k = omega/c*ones(1,N);

for n = 1:N
    sin_theta = sqrt(eps(1)*mu(1))*sin(theta_inc)/sqrt(eps(n)*mu(n));
    cos_theta(n) = sqrt( 1 - sin_theta^2 );
    if( sin_theta > 1 )
        k(n) = 1i*k(n);
    end
end

% расчет полного поля
f = zeros(4*N,1);
M = zeros( 4*N );

M(1,2) = 1i*sqrt(eps(1)*mu(1));
M(1,4) = 1i*mu(1);

for n = 1:N
    M(4*n-2,4*n-3) =  1i*k(n)*c*zmid(n);
    M(4*n-2,4*n-2) =  1i*k(n)*c;
    M(4*n-2,4*n-1) = -c/eps(n);
    
    M( 4*n-1 ,4*n-3) =  c/mu(n);
    M( 4*n-1 ,4*n-1) = -1i*k(n)*c*zmid(n);
    M( 4*n-1 , 4*n ) = -1i*k(n)*c;
end
if( polarization == 1 )
    for n = 2:N    
        M(4*n-4,4*n-7) = -z(n);
        M(4*n-4,4*n-6) = -1;
        M(4*n-4,4*n-3) =  z(n);
        M(4*n-4,4*n-2) =  1;
        
        M(4*n-3,4*n-5) =  z(n)*cos_theta(n-1);
        M(4*n-3,4*n-4) =    1 *cos_theta(n-1);
        M(4*n-3,4*n-1) = -z(n)*cos_theta( n );
        M(4*n-3, 4*n ) = -  1 *cos_theta( n );
    end
else
    for n = 2:N
        M(4*n-4,4*n-7) = -z(n)*cos_theta(n-1);
        M(4*n-4,4*n-6) = -  1 *cos_theta(n-1);
        M(4*n-4,4*n-3) =  z(n)*cos_theta( n );
        M(4*n-4,4*n-2) =    1 *cos_theta( n );
        
        M(4*n-3,4*n-5) =  z(n);
        M(4*n-3,4*n-4) =  1;
        M(4*n-3,4*n-1) = -z(n);
        M(4*n-3, 4*n ) = -1;
    end
end

M(4*N,4*N-3) = -1i*sqrt(eps(N)*mu(N))*a;
M(4*N,4*N-2) = -1i*sqrt(eps(N)*mu(N));
M(4*N,4*N-1) =  1i*mu(N)*a;
M(4*N, 4*N ) =  1i*mu(N); 

f(1) = 2*1i*sqrt(eps(1)*mu(1))*E0;
for n = 1:N
    f(4*n-2) =  4*pi*j(n)/eps(n);    
end
for n = 2:N
    f(4*n-3) = -4*pi*j_surf(n);
end
f(4*N) = 0;

% M = zeros( 2*N );
% f = zeros(2*N,1);
% 
% M(1,1) = -1i*k(1)*zmid(1)^2*mu(1) + 1/1i/k(1)/eps(1);
% M(1,2) = -1i*k(1)*zmid(1)*mu(1) + sqrt( mu(1)/eps(1) );
% f(1)   = 2*E0 - 4*pi*j(1)/c/eps(1)/1i/k(1);
% 
% for n = 2:N
%     M(2*n-2,2*n-3) = 1/k(n-1)/eps(n-1) - ( z(n)-zmid(n-1) )*k(n-1)*mu(n-1)*zmid(n-1);
%     M(2*n-2,2*n-2) = k(n-1)*mu(n-1)*( z(n)-zmid(n-1) );
%     M(2*n-2,2*n-1) = 1/k( n )/eps( n ) - ( z(n)-zmid( n ) )*k( n )*mu( n )*zmid( n );
%     M(2*n-2, 2*n ) = k( n )*mu( n )*( z(n)-zmid( n ) );
%     f(2*n-2) = 4*pi/c*j(n-1)/eps(n-1)/k(n-1) - 4*pi/c*j(n)/eps(n)/k(n);
%     
%     M(2*n-1,2*n-3) =  z(n)*cos_theta(n-1);
%     M(2*n-1,2*n-2) =  cos_theta(n-1);
%     M(2*n-1,2*n-1) = -z(n)*cos_theta( n );
%     M(2*n-1, 2*n ) = -cos_theta( n );
%     f(2*n-1) = 0;%4*pi*j_surf(n)/c;
% end
% 
% M(2*N,2*N-1) = 1/1i/k(N)/eps(N) + ( a-zmid(N) )*1i*k(N)*mu(N)*zmid(N) - a*sqrt( mu(N)/eps(N) );
% M(2*N, 2*N ) = (a - zmid(N)) * 1i*k(N)*mu(N) - sqrt(mu(N)/eps(N));
% f(2*N) = -4*pi*j(N)/c/eps(N)/1i/k(N);

x = M\f;

Fields = zeros(2,N);
for n = 1:N
    Fields(1,n) = x(4*n-3)*z(n) + x(4*n-2);
    Fields(2,n) = x(4*n-1)*z(n) + x( 4*n );
end

% A = zeros(1,N);
% B = zeros(1,N);
% C = zeros(1,N);
% D = zeros(1,N);
% for n = 1:N
%     C(n) = x(2*n-1);
%     D(n) = x( 2*n );
%     A(n) = 1i*k(n)*mu(n)*( C(n)*zmid(n) + D(n) );
%     B(n) = C(n)/1i/k(n)/eps(n) - A(n)*zmid(n) + 4*pi*j(n)/1i/eps(n)/k(n)/c;
% end
% for n = 1:N
%     Fields(1,n) = A(n)*z(n) + B(n);
%     Fields(2,n) = C(n)*z(n) + D(n);
% end

Fields_before = zeros(2,N1);
for n = 1:N1
    Fields_before(:,n) = Fields(:,n);
end

% расчет падающего поля
f = zeros(4*N1,1);
M = zeros( 4*N1 );

M(1,2) = 1i*sqrt(eps(1)*mu(1));
M(1,4) = 1i*mu(1);

for n = 1:N1
    M(4*n-2,4*n-3) =  1i*k(n)*c*zmid(n);
    M(4*n-2,4*n-2) =  1i*k(n)*c;
    M(4*n-2,4*n-1) = -c/eps(n);
    
    M( 4*n-1 ,4*n-3) =  c/mu(n);
    M( 4*n-1 ,4*n-1) = -1i*k(n)*c*zmid(n);
    M( 4*n-1 , 4*n ) = -1i*k(n)*c;
end
if( polarization == 1 )
    for n = 2:N1
        M(4*n-4,4*n-7) = -z(n);
        M(4*n-4,4*n-6) = -1;
        M(4*n-4,4*n-3) =  z(n);
        M(4*n-4,4*n-2) =  1;
        
        M(4*n-3,4*n-5) =  z(n)*cos_theta(n-1);
        M(4*n-3,4*n-4) =    1 *cos_theta(n-1);
        M(4*n-3,4*n-1) = -z(n)*cos_theta( n );
        M(4*n-3, 4*n ) = -  1 *cos_theta( n );
    end    
else
    for n = 2:N1
        M(4*n-4,4*n-7) = -z(n)*cos_theta(n-1);
        M(4*n-4,4*n-6) = -  1 *cos_theta(n-1);
        M(4*n-4,4*n-3) =  z(n)*cos_theta( n );
        M(4*n-4,4*n-2) =    1 *cos_theta( n );
        
        M(4*n-3,4*n-5) =  z(n);
        M(4*n-3,4*n-4) =  1;
        M(4*n-3,4*n-1) = -z(n);
        M(4*n-3, 4*n ) = -1;
    end
end

M(4*N1,4*N1-3) = -1i*sqrt(eps(1)*mu(1))*a;
M(4*N1,4*N1-2) = -1i*sqrt(eps(1)*mu(1));
M(4*N1,4*N1-1) =  1i*mu(1)*a;
M(4*N1, 4*N1 ) =  1i*mu(1); 

f(1) = 2*1i*sqrt(eps(1)*mu(1))*E0;
for n = 1:N1
    f(4*n-2) =  4*pi*j(n)/eps(n);    
end
for n = 2:N1
    f(4*n-3) = -4*pi*j_surf(n);
end
f(4*N1) = 0;

x1 = M\f;

% M = zeros( 2*N1 );
% f = zeros(2*N1,1);
% 
% M(1,1) = -1i*k(1)*zmid(1)^2*mu(1) + 1/1i/k(1)/eps(1);
% M(1,2) = -1i*k(1)*zmid(1)*mu(1) + sqrt( mu(1)/eps(1) );
% f(1)   = 2*E0 - 4*pi*j(1)/c/eps(1)/1i/k(1);
% 
% for n = 2:N1
%     M(2*n-2,2*n-3) = 1/k(n-1)/eps(n-1) - ( z(n)-zmid(n-1) )*k(n-1)*mu(n-1)*zmid(n-1);
%     M(2*n-2,2*n-2) = k(n-1)*mu(n-1)*( z(n)-zmid(n-1) );
%     M(2*n-2,2*n-1) = 1/k( n )/eps( n ) - ( z(n)-zmid( n ) )*k( n )*mu( n )*zmid( n );
%     M(2*n-2, 2*n ) = k( n )*mu( n )*( z(n)-zmid( n ) );
%     f(2*n-2) = 4*pi/c*j(n-1)/eps(n-1)/k(n-1) - 4*pi/c*j(n)/eps(n)/k(n);
%     
%     M(2*n-1,2*n-3) =  z(n)*cos_theta(n-1);
%     M(2*n-1,2*n-2) =  cos_theta(n-1);
%     M(2*n-1,2*n-1) = -z(n)*cos_theta( n );
%     M(2*n-1, 2*n ) = -cos_theta( n );
%     f(2*n-1) = 0;%4*pi*j_surf(n)/c;
% end
% 
% M(2*N1,2*N1-1) = 1/1i/k(N1)/eps(N1) + ( a-zmid(N1) )*1i*k(N1)*mu(N1)*zmid(N1) - a*sqrt( mu(N1)/eps(N1) );
% M(2*N1, 2*N1 ) = (a - zmid(N1)) * 1i*k(N1)*mu(N1) - sqrt(mu(N1)/eps(N1));
% f(2*N1) = -4*pi*j(N1)/c/eps(N1)/1i/k(N1);
% 
% x1 = M\f;


Fields_inc = zeros(1,N1);
for n = 1:N1
    Fields_inc(1,n) = x1(4*n-3)*z(n) + x1(4*n-2);
    Fields_inc(2,n) = x1(4*n-1)*z(n) + x1( 4*n );
end

% A = zeros(1,N);
% B = zeros(1,N);
% C = zeros(1,N);
% D = zeros(1,N);
% for n = 1:N1
%     C(n) = x1(2*n-1);
%     D(n) = x1( 2*n );
%     A(n) = 1i*k(n)*mu(n)*( C(n)*zmid(n) + D(n) );
%     B(n) = C(n)/1i/k(n)/eps(n) - A(n)*zmid(n) + 4*pi*j(n)/1i/eps(n)/k(n)/c;
% end
% for n = 1:N1
%     Fields_inc(1,n) = A(n)*z(n) + B(n);
%     Fields_inc(2,n) = C(n)*z(n) + D(n);
% end

% расчет отраженного поля
Fields_refl = Fields_before - Fields_inc;

I_refl = Fields_refl(1,:).*conj(Fields_refl(1,:));
I_refl_norm = I_refl/abs(E0)^2;
reflection = sqrt( sum(I_refl.^2)/length(I_refl) );
reflection_norm = sqrt( sum(I_refl_norm.^2)/length(I_refl_norm) );

% расчет прошедшего поля
% E_tr = zeros(1,N2);
% H_tr = zeros(1,N2);
% for n = 1:N2
%     E_tr(n) = E(n+N-N2);
%     H_tr(n) = H(n+N-N2);
% end
% I_tr = abs(E_tr).^2/abs(E0)^2;
% transmittance = sqrt( sum( I_tr.^2 )/length(I_tr) );
transmittance = abs(Fields(1,end))^2/abs(E0)^2;
end