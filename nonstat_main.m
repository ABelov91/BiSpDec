% [l] = 1 мкм = 1е-6 м
% [t] = 10 фс = 1е-14 с
% [v] = 1е8 м/с = 1/3 c0
% [omega] = 1e+14 с^{-1}

global theta_inc;
global omega_mid;
global polarization;
global fluct;
global problem;
global thick1;
global thick2;

tic;

c0 = 2.997924580; % скорость импульса в вакууме

problem = 2; % 0 -> тестовая, 1 -> старая (без Ag), 2 -> ФК 1 (160/112)*7+260, 3 -> ФК 2 (207/145)*7+185, 4 -> задание толщин "снаружи"
% 5 -> ФК1 в обратную сторону 260+(112/160)*7
theta_degrees = 45;
theta_radians = pi*(theta_degrees/180);
theta_inc = asin( 1.51/1.49*sin( theta_radians ) ); % угол падения от нормали
% theta_inc = theta_radians; % угол падения от нормали
polarization = 1; % 1 -> s, 2 -> p
step = 1; % шаг по частоте: 0 -> равномерный, 1 -> автоматический

% падающий импульс (гауссов)
lambda0 = 0.800; % центральная длина волны
delta_lambda_pulse = 0.032; % ширина
lambda_min = lambda0 - 10*delta_lambda_pulse;
lambda_max = lambda0 + 10*delta_lambda_pulse;

omega_max = 2*pi*c0/lambda_min;
omega_min = 2*pi*c0/lambda_max;
omega_mid = 2*pi*c0/lambda0;

T1 = -5;
T2 =  85;

materials;

pair_num = 7;
fluct = zeros(1,2*pair_num+1);
[z,h,~,Nums] = geometry;

N0 = length(h);
P0 = 100;
Q0 = 100;
h0 = ( omega_max - omega_min )/50;
thick_num = 2; % число сеток
est_E     = zeros(1, thick_num-1);
est_E_r   = zeros(1, thick_num-1);
est_CF    = zeros(1, thick_num-1);
est_CF_C  = zeros(1, thick_num-1);
R_min_arg = zeros(1, thick_num);
R_min_val = zeros(1, thick_num);
CF_mesh   = zeros(thick_num-1,P0*2^(thick_num-1)+1);
t2_mesh   = zeros(thick_num-1,P0*2^(thick_num-1)+1);
weights   = zeros(thick_num-1,P0*2^(thick_num-1)+1);
nodes     = zeros(1,thick_num);
nodes_rich = zeros(1,thick_num-1);

for mesh = 1:thick_num
    P = P0*2^(mesh-1); % время
    Q = Q0*2^(mesh-1); % частота
    N = N0*2^(mesh-1); % пространство
    
    nodes(mesh) = (P*Q*N)^(1/3);
    
    tau = (T2-T1)/P;
    t  = zeros(1,P);
    for p = 1:P+1
        t(p)  =  T1 + tau*(p-1);
    end
    
    F = zeros(1,Q+1);
    
    if( mesh > 1 )
        z = zeros(1,N+1);
        h = zeros(1, N );

        for n = 1:N/2
            h(2*n-1) = h_prev(n)/2;
            h( 2*n ) = h_prev(n)/2;
        end

        for n = 1:N/2
            z(2*n-1) = z_prev(n);
            z( 2*n ) = 0.5*z_prev(n) + 0.5*z_prev(n+1);
        end
        z(N+1) = z_prev(N/2+1);
    end

    % решение набора стационарных задач и обратное преобразование
    E = zeros(P+1,N);
    H = zeros(P+1,N);
    
    N1 = Nums(1)*2^(mesh-1);
    N2 = Nums(6)*2^(mesh-1);
    
    E_r = zeros(P+1,N1);
    H_r = zeros(P+1,N1);
    
    reflection      = zeros(1,Q+1);
    reflection_norm = zeros(1,Q+1);
    transmittance   = zeros(1,Q+1);
    omega   = zeros(1,Q+1);
    h_omega = zeros(1,Q);

    h0 = h0/2;
    h_omega(1) = h0;
    omega_c = omega_min;
    
    omega(1) = omega_c;
    
    lam = 2*pi*c0/omega(1);
    F(1) = exp( -(lam-lambda0).^2/delta_lambda_pulse^2 );
    param = [N1, N2, F(1), 0, omega_c];
    [Flds_p, Flds_r_p, R, R_norm, T, cos_theta] = stationary_1d_single_mesh(param,z);
    reflection(1) = R;
    reflection_norm(1) = R_norm;
    transmittance(1) = T;
    
    if( step == 1 )
        h_omega(1) = 1/( 1 + 1000*abs( abs( 1 - R / (abs( F(1) ))^2 ) ) )*h0;
    else
        h_omega(1) = h0;
    end
    
    omega_p  = omega_c;
    omega_c  = omega_c + h_omega(1);
    omega(2) = omega_c;
  
    for q = 2:Q+1
        
        lam = 2*pi*c0/omega(q);
        F(q) = exp( -(lam-lambda0).^2/delta_lambda_pulse^2 );
        
        param = [N1, N2, F(q), 0, omega_c];
        [Flds, Flds_r, R, R_norm, T, cos_theta] = stationary_1d_single_mesh(param,z);
        
        reflection(q) = R;
        reflection_norm(q) = R_norm;
        transmittance(q) = T;

        h_omega_c = h_omega(q-1);
        for p = 1:P+1
            E(p,:)   = E(p,:)   + 0.5*( Flds(1,:) * exp( -1i*omega_c*t(p) ) + Flds_p(1,:) * exp( -1i*omega_p*t(p) ) )*h_omega_c/sqrt(2*pi);
            E_r(p,:) = E_r(p,:) + 0.5*( Flds_r(1,:)*exp( -1i*omega_c*t(p) ) + Flds_r_p(1,:)*exp( -1i*omega_p*t(p) ) )*h_omega_c/sqrt(2*pi);
        end
        
        if( step == 1 )
            h_omega(q) = 1/( 1 + 1000*abs( 1 - R_norm ) )*h0;
        else
            h_omega(q) = h0;
        end
        
        Flds_p   = Flds;
        Flds_r_p = Flds_r;
        
        omega_p = omega_c;
        
        omega_c = omega_c + h_omega(q);
        if( q <= Q )
            omega(q+1) = omega_c;
        end
    end
    [R_min_val(mesh),n0] = min(reflection_norm);
    R_min_arg(mesh) = 2*pi*c0/omega(n0);
    
    t1 = zeros(1,3*P+1);
    for q1 = 1:3*P+1
        t1(q1) = 2*T1 - T2 + (q1-1)*tau;
    end
    E_inc = zeros(1,3*P+1);
    for q2 = 1:Q+1
        if( q2 > 1 )
            h_omega_c = h_omega(q2-1);
            for q1 = 1:3*P+1
                E_inc(q1) = E_inc(q1) + ...
                    ( F( q2 )*exp( -1i*omega( q2 )*t1(q1) ) + ...
                      F(q2-1)*exp( -1i*omega(q2-1)*t1(q1) ) )*h_omega_c/sqrt(2*pi);
            end
        end
    end
    
    CF = zeros(1,2*P+1);   
    for q1 = 1:2*P+1
        for q2 = 1:P
            tmp1 = abs(E_r( q2 ,1))^2*abs( E_inc(q1+q2-1) )^2;
            tmp2 = abs(E_r(q2+1,1))^2*abs( E_inc( q1+q2 ) )^2;
            CF(q1) = CF(q1) + 0.5*(tmp1+tmp2)*tau;
        end
    end
    CF = CF/max(CF);
    
    if( mesh > 1 )
        E_cut   = zeros(P/2+1,N/2);
        E_r_cut = zeros(P/2+1,N1/2);
        CF_cut  = zeros(1,P+1);
        for n = 1:N/2
            for p = 1:P/2+1
                E_cut(p,n) = E(2*p-1,2*n-1);
            end
        end
        for n = 1:N1/2
            for p = 1:P/2+1
                E_r_cut(p,n) = E_r(2*p-1,2*n-1);
            end
        end
        for p = 1:P+1
            CF_cut(p) = CF(2*p-1);
        end
        
        num = (2^thick_num - 2^mesh)*P0/2;
        tmp = [CF_cut,  zeros(1,num) ];
        CF_mesh(mesh-1,:) = tmp;
        
        for q1 = 1:P+1
            t2_mesh(mesh-1,q1) = (T1-T2) + (q1-1)*2*tau;
        end
        
        rich_E   = ( E_cut   - E_prev )  /3;        
        rich_E_r = ( E_r_cut - E_r_prev )/3;
        rich_CF  = ( CF_cut - CF_prev ) /3;
        
        est_E(mesh-1)    = sqrt( sum( abs( rich_E ).^2, 'all' )/(0.5*N)/(0.5*P+1) );
        est_E_r(mesh-1)  = sqrt( sum( abs(rich_E_r).^2, 'all' )/(0.5*N)/(0.5*P+1) );
        est_CF(mesh-1)   = sqrt( rich_CF*rich_CF'/(P+1) );
        est_CF_C(mesh-1) = max( abs( rich_CF ) );
        
        nodes_rich(mesh-1) = nodes(mesh);
        
        tmp = [rich_CF,  zeros(1,num) ];
        weights(mesh-1,:) = abs( tmp );
    end
    
    E_prev   = E;
    E_r_prev = E_r;
    z_prev   = z;
    h_prev   = h;
    CF_prev  = CF;
    
    lambda = 2*pi*c0./omega;
end

extrapolation;
illustrations;
output;
file_name = num2str(theta_degrees);
file_name = strcat(file_name,'.mat');
save(file_name,'-v7.3');

toc