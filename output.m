angle = num2str(theta_degrees);
wavelength = num2str(lambda0);
if(problem == 2)
    file_name = strcat('PC1_',angle);
elseif( problem == 3 )
    file_name = strcat('PC2_',angle);
else
    file_name = strcat('Other_',angle);
end
% file_name = strcat(file_name,'_');
% file_name = strcat(file_name,wavelength);
file_name = strcat(file_name,'_output.txt');
fid0 = fopen(file_name,'w');
fprintf(fid0,'Thickening meshes output\r\n');
fprintf(fid0,'lambda0 = %.3f, delta_lambda = %.3f\r\n\n', lambda0, delta_lambda_pulse);
fprintf(fid0,'LIFETIMES:\r\n');
fprintf(fid0,'%.4f fs\t', lifetime);
fprintf(fid0,'\r\n\n');

fprintf(fid0,'Error estimations:\r\n');
fprintf(fid0,'%.4f fs,\t', est_lifetime);
fprintf(fid0,'\r\n\n');

fprintf(fid0,'Extrapolated value:\r\n');
fprintf(fid0,'%.4f fs,\t', extrap_lifetime);
fprintf(fid0,'\r\n\n');

fprintf(fid0,'Mean deviation of the approximation:\r\n');
fprintf(fid0,'%.4e\t', r1);
fprintf(fid0,'\r\n\n');

fprintf(fid0,'Mean-square deviation of the approximation:\r\n');
fprintf(fid0,'%.4e\t', r2);
fprintf(fid0,'\r\n\n');

fprintf(fid0,'LG ERRORS\r\n');
fprintf(fid0,'Cross-corr:\r\n');
fprintf(fid0,'%.4f\t', log10(est_CF));
fprintf(fid0,'\r\n\n');

fprintf(fid0,'Total field:\r\n');
fprintf(fid0,'%.4f\t', log10(est_E));
fprintf(fid0,'\r\n\n');

fprintf(fid0,'Reflected field:\r\n');
fprintf(fid0,'%.4f\t', log10(est_E_r));
fprintf(fid0,'\r\n\n');

fprintf(fid0,'REFLECTION SPECTRUM\r\n');
fprintf(fid0,'resonance position:\r\n');
fprintf(fid0,'%.6f\t', R_min_arg);
fprintf(fid0,'\r\n\n');

fprintf(fid0,'resonance depth:\r\n');
fprintf(fid0,'%.6f\t', R_min_val);
fprintf(fid0,'\r\n\n');

fclose(fid0);
