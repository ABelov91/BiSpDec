if( thick_num > 1 )
    t_bsw_min = 10;
    t_bsw_max = 80;

    lifetime = zeros(1,thick_num-1);
    tail_intens = zeros(1,thick_num-1);
    est_lifetime    = zeros(1,thick_num-2);
    extrap_lifetime = zeros(1,thick_num-2);
    delta_lifetime  = zeros(1,thick_num-1);
    est_tail_intens    = zeros(1,thick_num-2);
    extrap_tail_intens = zeros(1,thick_num-2);
    est_rmin_val    = zeros(1,thick_num-2);
    extrap_rmin_val = zeros(1,thick_num-2);
    est_rmin_arg    = zeros(1,thick_num-2);
    extrap_rmin_arg = zeros(1,thick_num-2);
    r1 = zeros(1,thick_num-1);
    r2 = zeros(1,thick_num-1);

    figure;
    for mesh = 1:thick_num-1
        xi = double.empty;
        eta = double.empty;
        rho = double.empty;

        count = 1;

        for n = 1:P0*2^(mesh-1)
            if( t_bsw_min <= -t2_mesh(mesh,n) && -t2_mesh(mesh,n) <= t_bsw_max )
                xi(count) = -t2_mesh(mesh,n);
                eta(count) = log10(CF_mesh(mesh,n));
                rho(count) = 1/weights(mesh,n)^2;
                count = count + 1;
            end
        end
        points_num = count - 1;

        Mxi  = sum(xi.*rho)/sum(rho);
        
        A = sum(rho.*eta)/sum(rho);
        delta_A = 1/sum(rho);
        
        P1 = xi - Mxi*ones(1,points_num);
        B = sum( rho.*eta.*P1 )/sum( P1.^2.*rho );
        delta_B = 1/sum( P1.^2.*rho );

        r1(mesh) = sum( eta - A - B*(xi - Mxi) )/points_num;
        r2(mesh) = sqrt( sum( (eta - A - B*(xi - Mxi)).^2 )/points_num );

        hold on;
        plot(-t2_mesh(mesh,:),log10(CF_mesh(mesh,:)))
        plot(xi,A + B*(xi - Mxi),'-g')
        
        tail_intens(mesh) = max( A + B*(xi - Mxi) );
        lifetime(mesh) = -10/B/log(10);
        delta_lifetime(mesh) = delta_B/B^2*10/log(10);
        
        if(mesh > 1)
            est_lifetime(mesh-1)    = ( lifetime(mesh) - lifetime(mesh-1) )/3;
            extrap_lifetime(mesh-1) = lifetime(mesh) + est_lifetime(mesh-1);
            est_tail_intens(mesh-1)    = ( tail_intens(mesh) - tail_intens(mesh-1) )/3;
            extrap_tail_intens(mesh-1) = tail_intens(mesh) + est_tail_intens(mesh-1);
            est_rmin_val(mesh-1)    = ( R_min_val(mesh) - R_min_val(mesh-1) )/3;
            extrap_rmin_val(mesh-1) = R_min_val(mesh) + est_rmin_val(mesh-1);
            est_rmin_arg(mesh-1)    = ( R_min_arg(mesh) - R_min_arg(mesh-1) )/3;
            extrap_rmin_arg(mesh-1) = R_min_arg(mesh) + est_rmin_arg(mesh-1);
        end
    end
%     xlabel('t \times 10^{-1}, fs')
%     ylabel('lg CF, arb. units')
end