% https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0


% Define the shape and scale parameters of the gamma distribution
shape = 2;
scale = 3;

% Generate random data from the gamma distribution
data = gamrnd(shape, scale, [1, 1000]);

% % Plot histogram of the generated data
% hist(data);
% % Plot the PDF of the gamma distribution with the same shape and scale parameters
% x = 0:0.1:20;
% pdf = gampdf(x, shape, scale)*1000;
% hold on
% plot(x, pdf, 'r');
% legend('Simulated Data', 'Gamma PDF');

% Define parameters for negative binomial distribution
r = 10;
p = 0.1;
samples = nbinrnd(r,p,1,1000);




