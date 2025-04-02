clear; clc; close all;

syms omega k nu sigma;

rho = 1e3; % kg/m^3
g = 9.81; % m/s^2

% full dispersion equation
F(omega, k, nu, sigma) = ((-1i*omega+2*nu*k^2)^2+g*k+sigma*k^3/rho)^2-(4*nu^2*k^3)^2*(k^2-1i*omega/nu); % sq terms
% F(omega, k, nu, sigma) = (((-1i*omega+2*nu*k^2)^2+g*k+sigma*k^3/rho)-(4*nu^2*k^3)*sqrt(k^2-1i*omega/nu)); % original
% disp('dispersion equation'); disp(F);

% group velocity
dw_dk = -diff(F, k)/diff(F, omega); % using implicit function derivative property

% f_sweep_list = [10, 20, 30, 40, 50, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 40, 40, 40, 40, 40]; % Hz
% sigma_sweep_list = [70, 70, 70, 70, 70, 20, 70, 100, 200, 500, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70]*1e-3; % N/m
% mu_sweep_list = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 10, 25, 50, 100, 1, 10, 25, 50, 100]*1e-3; % Pa.s


f_sweep_list = [20];
sigma_sweep_list = [70]*1e-3;
mu_sweep_list = [150]*1e-3;


for id=1:length(f_sweep_list)
    
    f_val = f_sweep_list(id);
    sigma_val = sigma_sweep_list(id);
    mu_val = mu_sweep_list(id);
    omega_val = 2*pi*f_val; % rad/s
    nu_val = mu_val/rho; % SI unit


    % [coefficients, powers] = coeffs(F(omega_val, k, nu_val, sigma_val), k);
    % for i=1:length(coefficients)
    %     disp(['power of k = ', string(powers(i)), 'coeff of k = ', num2str(double(coefficients(i)), "%0.1d")]);
    % end

    disp(['f = ', num2str(f_val), ' Hz, sigma = ', num2str(sigma_val*1e3), ' mN/m, mu = ', num2str(mu_val*1e3), ' mPa.s']);

    solutions = vpasolve(F(omega_val, k, nu_val, sigma_val)==0);
    % disp('all solutions of k = '); disp(double(solutions));
    disp(['total number of solutions = ', num2str(length(solutions))]);
    disp('valid solutions');

    k_list = nan(length(solutions), 1);
    sigma_list = nan(length(solutions), 1);
    mu_calc1_list = nan(length(solutions), 1);
    mu_calc2_list = nan(length(solutions), 1);

    for i = 1:length(solutions)
        m = sqrt(solutions(i)^2 - 1i*omega_val/nu_val);
        % positive value of attenuation coeff is valid as it results in decay because of e^(ikx-wt) term

        if(real(solutions(i))>0) && (imag(solutions(i))>0) && (real(m)>0)

            %            disp(solutions(i));
               k_val = double(solutions(i));

               lambda = double(2*vpa(pi)/real(solutions(i))); % m
               error = abs(double((subs(F(omega_val, k_val, nu_val, sigma_val))))); % abs() because of complex number

               phase_vel1 = real(omega_val/k_val);% % m/s
               phase_vel2 = double(f_val*lambda); % m/s
               group_vel = real(double((subs(dw_dk(omega_val, k_val, nu_val, sigma_val))))); % m/s, real part of dw/dk evaluated at calculated k and omega

               sigma_calc = double((omega_val^2/real(k_val)-g)*rho/real(k_val)^2); % using Eq. 2 from paper, dispersion equation for inviscid liquid
               mu_calc1 = real(3*rho*phase_vel1*imag(k_val)/4/real(k_val)^2); % using Eq. 3 from paper
               % solving actual dispersion relation (sq root terms) to get single mu values as output
    %            mu_calc2 = real(double(vpasolve(sqrt(F(omega_val, k_val, nu, sigma_calc))==0)*rho));
               mu_calc2 = real(double(vpasolve(F(omega_val, k_val, nu, sigma_calc)==0))*rho);

               output1 = ['lambda = ', num2str(lambda*1e3, "%0.1f"), ' mm,'];
               output2 = [' k = ', num2str(k_val, "%0.1f"), ' 1/m,']; 
               output3 = [' error = ', num2str(error, "%0.1d")];

               output4 = ['phase vel (w/k) = ', num2str(phase_vel1*1e3, "%0.1f"), ' mm/s,'];
               output5 = [' phase vel (f*lambda) = ', num2str(phase_vel2*1e3, "%0.1f"), ' mm/s,'];
               output6 = [' group speed (dw/dk) = ', num2str(group_vel*1e3, "%0.1f"), ' mm/s,'];

               output7 = ['sigma calc = ', num2str(sigma_calc*1e3, "%0.1f"), ' mN/m,'];
               output8 = [' mu calc (simple eq.) = ', num2str(mu_calc1*1e3, "%0.1f"), ' mPa.s,'];
               output9 = [' mu calc (full disp eq.) = ', num2str(mu_calc2'*1e3, '%0.1f,'), ' mPa.s'];

               disp([output1, output2, output3]);
               disp([output4, output5, output6]);
               disp([output7, output8, output9]);
    %            disp(['mu calc2 = ', num2str(mu_calc2'*1e3, '%0.1f,'), ' mPa.s']);


               k_list(i) = k_val;
               sigma_list(i) = sigma_calc*1e3;
               mu_calc1_list(i) = mu_calc1*1e3;
    %            mu_calc2_list(i) = mu_calc2*1e3;
               disp('--------');         
        end
    end

    % disp(k_list');
    % disp(sigma_list');
    % disp(mu_calc1_list');
    % disp(mu_calc2_list');

    % disp(['mean k=', num2str(mean(k_list, 'omitnan'), '%0.1f'), ', mean sigma=', num2str(mean(sigma_list, 'omitnan'), '%0.1f'), ', mean mu_calc2=', num2str(mean(mu_calc2_list, 'omitnan'), '%0.1f'), ', mean mu_calc1=', num2str(mean(mu_calc1_list, 'omitnan'), '%0.1f')]);
    disp(['mean k=', num2str(mean(k_list, 'omitnan'), '%0.1f'), ', mean sigma=', num2str(mean(sigma_list, 'omitnan'), '%0.1f'), ', mean mu_calc1=', num2str(mean(mu_calc1_list, 'omitnan'), '%0.1f')]);
    disp('------------------------------------------------------');
end
%-------------------
% simplified dispersion eq, behroozi 2011 paper
% c0 = -omega^2;
% c1 = g;
% c2 = 0;
% c3 = sigma/rho-sqrt(8*nu^3*omega);
% c4 = 4*nu^2;
% poly_coeff = [c4 c3 c3 c2 c1 c0];
% disp(newline); disp('poly_coeff = '); disp(poly_coeff);j*
% sol = roots(poly_coeff);
% disp(newline); disp('simplified solution (behroozi 2011), k = '); 
% disp(sol);