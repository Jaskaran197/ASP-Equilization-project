

close all
clc
rng('shuffle')

N = 600;
K = 300;

% effect of eigenvalue spread
h1 = [0.2194, 1.0, 0.2194];
M=9;
sigma_noise=0.01;

step_size = 0.075;
e = run_ensemble(N,K,h1,sigma_noise,M,step_size,false);
figure('DefaultAxesFontSize',24);
semilogy(e,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");
hold on

h2=[0.2798,1.0,0.2194];
e = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);
semilogy(e,'LineWidth',2);

h2=[0.3365,1.0,0.3365];
e = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);
semilogy(e,'LineWidth',2);
legend('Channel 1 - [0.2194, 1.0, 0.2194]', 'Channel 2 - [0.2798, 1.0, 0.2798]', 'Channel 3 - [0.3365, 1.0, 0.3365]')
title('Effect of Eigenvalue spread');
saveas(gcf,'Effect of eignevalue spread.png');
hold on


% effect of filter order
N = 800;
K = 400;
M = 9;
h2 = [0.2798, 1.0, 0.2798];
step_size = 0.075;
e = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);
figure('DefaultAxesFontSize',24);
semilogy(e,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");

hold on


M =11;
h2 = [0.2798, 1.0, 0.2798];
step_size = 0.075;
e = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);

semilogy(e,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");
hold on


M =21;
h2 = [0.2798, 1.0, 0.2798];
step_size = 0.075;
e = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);

semilogy(e,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");

hold on
legend('M-9','M-11','M-21')
title('Effect of filter order');
saveas(gcf,'Effect of filter order.png');

N = 800;
K = 400;
h1 =[0.2194, 1.0, 0.2194];
M=9;
sigma_noise=0.01;

% effect of filter order


h2 = [0.2798, 1.0, 0.2798];
step_size = 0.0375;
e = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);
figure('DefaultAxesFontSize',24);

semilogy(e,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");
hold on


M = 11;
h2 = [0.2798, 1.0, 0.2798];
step_size = 0.0375;
e = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);

semilogy(e,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");
hold on


M = 21;
h2 = [0.2798, 1.0, 0.2798];
step_size = 0.0375;
e = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);

semilogy(e,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");
legend('M-9','M-11','M-21')
title('Effect of filter order');
hold on





% effect of filter order

M = 9;
h2 = [0.2798, 1.0, 0.2798];
step_size = 0.0125;
e = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);

figure('DefaultAxesFontSize',24);
semilogy(e,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");
hold on


M = 11;
h2 = [0.2798, 1.0, 0.2798];
step_size = 0.0125;
e = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);

semilogy(e,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");
hold on


M = 21;
h2 = [0.2798, 1.0, 0.2798];
step_size = 0.0125;
e = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);

semilogy(e,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");
legend('M-9','M-11','M-21')
title('Effect of filter order');
hold on



channel1 = [0.2194, 1.0, 0.2194];
M=9;
step_size  = 0.0125;
N = 1600;
K = 300;

sigma_noise=0.01;



b = run_ensemble(N,K,channel1,sigma_noise,M,step_size,false);
figure('DefaultAxesFontSize',25);
semilogy(b,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");

hold on


channel1 = [0.2194, 1.0, 0.2194];
M=9;
step_size  = 0.025;
N = 1600;
K = 300;

sigma_noise=0.01;



bd = run_ensemble(N,K,channel1,sigma_noise,M,step_size,false);
%figure('DefaultAxesFontSize',25);
semilogy(bd,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");

hold on

channel1 = [0.2194, 1.0, 0.2194];
M=9;
step_size  = 0.075;
N = 1600;
K = 300;

sigma_noise=0.01;



bc = run_ensemble(N,K,channel1,sigma_noise,M,step_size,false);

semilogy(bc,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");
legend('0.0125','0.025','0.075')
title('Effect of step size');
saveas(gcf,'Effect of step size.png');
%figure('DefaultAxesFontSize',25);
hold on



%%% normalized mse vs standard

h2=[0.2798,1.0,0.2194];
M=9;
step_size  = 0.025;
N = 1600;
K = 300;

sigma_noise=0.01;


figure('DefaultAxesFontSize',25);
bc = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);
semilogy(bc,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");

hold on

step_size  = 0.075;
bc = run_ensemble(N,K,h2,sigma_noise,M,step_size,false);
semilogy(bc,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");

hold on
bc = run_ensemble(N,K,h2,sigma_noise,M,step_size,true);
semilogy(bc,'LineWidth',2);
xlabel('Time Sample (n)');
ylabel("Normalized MSE");


legend('Step size - 0.025','Step size - 0.075','Normalized MSE');
title('Standard LMS vs Normalized LMS');

saveas(gcf,'Standard LMS vs Normalized LMS.png');




function e = run_ensemble(N,K,h, sigma_noise, M, step, normalized_lms)
    e = zeros(K,N);
    for i = 1:K
        e(i,:) = run_experiment(N, h, sigma_noise, M, step, normalized_lms);
    end
    e = mean(e);

end





















function e = run_experiment(N, h, sigma_noise, M, step, normalized_lms)
    d = 2*((randn(1,N)>0)-0.5);
    x = conv(d,h);
    x = x(1:N);
    noise = sigma_noise*randn(1,N);
    u= x+noise;
    w = zeros(M,N);
    d_hat = zeros(1,N);
    e = zeros(1,N);

    delta1 = floor(0.5*length(h));
    delta2 = floor(0.5*M);
    delay = delta1  + delta2;

    for i = 1:N
        for j = 0:M-1
            if(i-j>0)
                d_hat(i) = d_hat(i) + w(j+1,i)*u(i-j);
            end
        end
        if (i>delay)
            e(i) = d(i-delay) - d_hat(i);
        end
        if(i<N)
            u_vector = zeros(M,1);
            u_vector(1:min(i,M)) = transpose(flip(u(max(1,i-M+1): i)));
            if(normalized_lms)
                step = 1/(norm(u_vector)).^2;
            end
            w(:,i+1) = w(:,i) + step*u_vector*e(i);
        end
    end
    e=e.^2;
end