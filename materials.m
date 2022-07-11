global ppSIO2;
global ppTA2O5;

materials_choice = 1; % 1 -> IZOVAC, 2 -> Rodriguez2016, 3 -> Rodriguez2016 & Bright2013, 4 -> Gao2012

c0 = 2.997924580;

switch materials_choice
    case 1
    mul_SIO2  = 1;
    mul_TA2O5 = 2.5;
    file_title = 'CRYSTALS_SiO2_IZOVAC.dat';
    fid = fopen(file_title,'r');
    S1 = fscanf(fid,'%g');
    fclose(fid);

    Num = length(S1)/3;
    lambda_SiO2 = zeros(1,Num);
    omega_SiO2 = zeros(1,Num);
    refr_index_SiO2 = zeros(1,Num);
    epsilon_SiO2 = zeros(1,Num);
    for n = 1:Num
        lambda_SiO2(n) = S1( 1 + (n-1)*3 )*1e-3;
        omega_SiO2(n) = 2*pi*c0/lambda_SiO2(n);
        refr_index_SiO2(n) = S1( 2 + (n-1)*3 ) + 1i*S1( 3 + (n-1)*3 )*mul_SIO2;
        epsilon_SiO2(n) = refr_index_SiO2(n)^2;
    end

    file_title = 'CRYSTALS_Ta2O5_IZOVAC.dat';
    fid = fopen(file_title,'r');
    S2 = fscanf(fid,'%g');
    fclose(fid);

    Num = length(S2)/3;
    lambda_Ta2O5 = zeros(1,Num);
    omega_Ta2O5 = zeros(1,Num);
    refr_index_Ta2O5 = zeros(1,Num);
    epsilon_Ta2O5 = zeros(1,Num);
    for n = 1:Num
        lambda_Ta2O5(n) = S2( 1 + (n-1)*3 )*1e-3;
        omega_Ta2O5(n) = 2*pi*c0/lambda_Ta2O5(n);
        refr_index_Ta2O5(n) = S2( 2 + (n-1)*3 ) + 1i*S2( 3 + (n-1)*3 )*mul_TA2O5;
        epsilon_Ta2O5(n) = refr_index_Ta2O5(n)^2;
    end
    
    case 2
    file_title = 'REFR_IND_n_SiO2_Rodriguez.txt';
    fid = fopen(file_title,'r');
    S1 = fscanf(fid,'%g');
    fclose(fid);

    file_title = 'REFR_IND_k_SiO2_Rodriguez.txt';
    fid = fopen(file_title,'r');
    S2 = fscanf(fid,'%g');
    fclose(fid);

    Num = length(S1)/2;
    lambda_SiO2 = zeros(1,Num);
    omega_SiO2 = zeros(1,Num);
    refr_index_SiO2 = zeros(1,Num);
    epsilon_SiO2 = zeros(1,Num);

    for n = 1:Num
        lambda_SiO2(n) = S1( 1 + (n-1)*2 );
        omega_SiO2(n) = 2*pi*c0/lambda_SiO2(n);
        refr_index_SiO2(n) = S1( 2 + (n-1)*2 ) + 1i*S2( 2 + (n-1)*2 );
        epsilon_SiO2(n) = refr_index_SiO2(n)^2;
    end

    file_title = 'REFR_IND_n_Ta2O5_Rodriguez.txt';
    fid = fopen(file_title,'r');
    S3 = fscanf(fid,'%g');
    fclose(fid);

    file_title = 'REFR_IND_k_Ta2O5_Rodriguez.txt';
    fid = fopen(file_title,'r');
    S4 = fscanf(fid,'%g');
    fclose(fid);

    Num = length(S3)/2;
    lambda_Ta2O5 = zeros(1,Num);
    omega_Ta2O5 = zeros(1,Num);
    refr_index_Ta2O5 = zeros(1,Num);
    epsilon_Ta2O5 = zeros(1,Num);

    for n = 1:Num
        lambda_Ta2O5(n) = S3( 1 + (n-1)*2 );
        omega_Ta2O5(n) = 2*pi*c0/lambda_Ta2O5(n);
        refr_index_Ta2O5(n) = S3( 2 + (n-1)*2 ) + 1i*S4( 2 + (n-1)*2 );
        epsilon_Ta2O5(n) = refr_index_Ta2O5(n)^2;
    end

    case 3
    file_title = 'REFR_IND_n_SiO2_Rodriguez.txt';
    fid = fopen(file_title,'r');
    S1 = fscanf(fid,'%g');
    fclose(fid);

    file_title = 'REFR_IND_k_SiO2_Rodriguez.txt';
    fid = fopen(file_title,'r');
    S2 = fscanf(fid,'%g');
    fclose(fid);

    Num = length(S1)/2;
    lambda_SiO2 = zeros(1,Num);
    omega_SiO2 = zeros(1,Num);
    refr_index_SiO2 = zeros(1,Num);
    epsilon_SiO2 = zeros(1,Num);

    for n = 1:Num
        lambda_SiO2(n) = S1( 1 + (n-1)*2 );
        omega_SiO2(n) = 2*pi*c0/lambda_SiO2(n);
        refr_index_SiO2(n) = S1( 2 + (n-1)*2 ) + 1i*S2( 2 + (n-1)*2 );
        epsilon_SiO2(n) = refr_index_SiO2(n)^2;
    end

    file_title = 'REFR_IND_n_Ta2O5_Bright.txt';
    fid = fopen(file_title,'r');
    S3 = fscanf(fid,'%g');
    fclose(fid);

    file_title = 'REFR_IND_k_Ta2O5_Bright.txt';
    fid = fopen(file_title,'r');
    S4 = fscanf(fid,'%g');
    fclose(fid);

    Num = length(S3)/2;
    lambda_Ta2O5 = zeros(1,Num);
    omega_Ta2O5 = zeros(1,Num);
    refr_index_Ta2O5 = zeros(1,Num);
    epsilon_Ta2O5 = zeros(1,Num);

    for n = 1:Num
        lambda_Ta2O5(n) = S3( 1 + (n-1)*2 );
        omega_Ta2O5(n) = 2*pi*c0/lambda_Ta2O5(n);
        refr_index_Ta2O5(n) = S3( 2 + (n-1)*2 ) + 1i*S4( 2 + (n-1)*2 );
        epsilon_Ta2O5(n) = refr_index_Ta2O5(n)^2;
    end

    case 4
    file_title = 'REFR_IND_n_SiO2_Gao.txt';
    fid = fopen(file_title,'r');
    S1 = fscanf(fid,'%g');
    fclose(fid);

    file_title = 'REFR_IND_k_SiO2_Gao.txt';
    fid = fopen(file_title,'r');
    S2 = fscanf(fid,'%g');
    fclose(fid);

    Num = length(S1)/2;
    lambda_SiO2 = zeros(1,Num);
    omega_SiO2 = zeros(1,Num);
    refr_index_SiO2 = zeros(1,Num);
    epsilon_SiO2 = zeros(1,Num);

    for n = 1:Num
        lambda_SiO2(n) = S1( 1 + (n-1)*2 );
        omega_SiO2(n) = 2*pi*c0/lambda_SiO2(n);
        refr_index_SiO2(n) = S1( 2 + (n-1)*2 ) + 1i*S2( 2 + (n-1)*2 );
        epsilon_SiO2(n) = refr_index_SiO2(n)^2;
    end

    file_title = 'REFR_IND_n_Ta2O5_Gao.txt';
    fid = fopen(file_title,'r');
    S3 = fscanf(fid,'%g');
    fclose(fid);

    file_title = 'REFR_IND_k_Ta2O5_Gao.txt';
    fid = fopen(file_title,'r');
    S4 = fscanf(fid,'%g');
    fclose(fid);

    Num = length(S3)/2;
    lambda_Ta2O5 = zeros(1,Num);
    omega_Ta2O5 = zeros(1,Num);
    refr_index_Ta2O5 = zeros(1,Num);
    epsilon_Ta2O5 = zeros(1,Num);

    for n = 1:Num
        lambda_Ta2O5(n) = S3( 1 + (n-1)*2 );
        omega_Ta2O5(n) = 2*pi*c0/lambda_Ta2O5(n);
        refr_index_Ta2O5(n) = S3( 2 + (n-1)*2 ) + 1i*S4( 2 + (n-1)*2 );
        epsilon_Ta2O5(n) = refr_index_Ta2O5(n)^2;
    end
end

ppSIO2 = pchip(omega_SiO2,epsilon_SiO2);
ppTA2O5 = pchip(omega_Ta2O5,epsilon_Ta2O5);