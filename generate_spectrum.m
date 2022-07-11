lambda0 = 0.8; % ‘ 1: 0.8268, 0.7297, 0.641294, 0.5565. ‘ 2: 0.8388, 0.750, 0.6848, 0.5566
delta_lambda = 0.005;
delta_lambda_pulse = 0.032;

lambda = 0.5:1e-3:1.2;

F = exp( -(lambda-lambda0).^2/delta_lambda_pulse^2 );
A = [lambda; F];

fid1 = fopen('spectrum_model.txt','w');
fprintf(fid1,'%3.8e %3.8e \r\n', A);
fclose(fid1);

figure; hold on;
plot(lambda,F)

lambda_max = max(lambda);
lambda_min = min(lambda);