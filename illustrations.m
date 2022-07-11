lambda = 2*pi*c0./omega;

figure;
hold on;
plot(lambda, reflection_norm,'-r')
% plot(lambda, reflection,'-r')
plot(lambda, transmittance,'-b')
% experiments;

% plot(lambda, log10(h_omega),'-or')
xlabel('\lambda, {\mu}m')
ylabel('Reflection')

% eps = zeros(1,N);
% mu  = zeros(1,N);
% z_left = zeros(1,N);
% z_mid = zeros(1,N);
% for n = 1:N
%     z_left(n) = z(n);
%     z_mid(n)  = 0.5*z(n) + 0.5*z(n+1);
%     eps(n) = diel(z_mid(n), omega_mid);
%     mu(n)  = magn(z_mid(n), omega_mid);
% end

% figure; %hold on;
% reduce = 5;
% P_temp = P/reduce;
% if(polarization == 1)
%     for p = 1:P_temp
%         count = ( p-1 )*reduce + 1;
%         plot(z_mid,real(eps),':k', z_mid,imag(eps),'--k',...
%              z_left,real(E_ans(count,:)),'-b', z_left,imag(E_ans(count,:)),'-g',...
%              z_left,real(H_ans(count,:).*cos_theta),'-r', z_left,imag(H_ans(count,:).*cos_theta),'-m')    
%         axis([0 z(end) -2 6])
%         pause(0.05)
%     end
% else
%     for p = 1:P_temp
%         count = ( p-1 )*reduce + 1;
%         plot(z_mid,real(eps),':k', z_mid,imag(eps),'--k',...
%              z_left,real(E_ans(count,:).*cos_theta),'-b', z_left,imag(E_ans(count,:).*cos_theta),'-g',...
%              z_left,real(H_ans(count,:)),'-r', z_left,imag(H_ans(count,:)),'-m')    
%         axis([0 z(end) -2 6])
%         pause(0.05)
%     end
% end

% figure; hold on;
% plot(log10(nodes_rich),log10(est_E),'-or')
% plot(log10(nodes_rich),log10(est_H),'-ob')
% legend( 'E Rich', 'H Rich' )
% xlabel('lg N')
% ylabel('lg error')

% figure;
% hold on;
% plot(lambda, transmittance,'-k')
% xlabel('\lambda, {\mu}m')
% ylabel('Transmittance')