function [x,mass,y0,cut_indices0,yfinal,cut_indices,I,iters] = runex(eg)
%RUNEX examples for computing local minimizers of the (MPPC) functional
%   Examples: (2,4,6,10 are examples from the paper)
%   1. Arc with noise in y-direction and straight line initial data
%   2. Spiral
%   3. Three spirals
%   4. Arc with spiral wrapped around it
%   5. Two concentric arcs
%   6. A noisy NxN grid with background noise
%   7. Cosmic data
%   8. Filamentary 3-d data
%   9. 'MPPC' letters with background noise
%   10. Zebrafish images


switch eg
    
    % Example 1: Arc with noise in y-direction and straight line initial data
    case 1
        i = (0:0.01:5)';
        n = length(i);
        x = [cos(i), sin(i)+0.5*rand(n,1)];
        j = -0.8:0.06:1;
        
        m = length(j);
        y0 = [j; 1.95*ones(1,m)]';
        cut_indices0 = [];
        mass = ones(1,n)./n;
        
        lambda1 = .002;
        rho = .4;
        tol = 10^-6;
        
        lambda2 = .05;
        max_m = 50;
        max_avg_turn = 100;
        normalize_data = 1;
        pause_bool = 0;
        
        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,x,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;
        
        
        % Example 2: Spiral
    case 2        
        n = 2000;
        df = @(s) sqrt(1+s^2);
        v = getUnifPoints(df,5*n,.001,3,14,.001);
        v = datasample(v,n);
        x_gen = [v.*cos(v),v.*sin(v)]; %generating curve
        x = x_gen + 1.5*randn(n,2);
        mass=ones(1,n)./n;
        
        %rng(1);
        y0 = [];
        cut_indices0 = [];
        
        lambda1 = .01;
        crit_dens = .07;
        lambda2 = 4/3*sqrt(lambda1/crit_dens);
        
        rho = .05;
        tol = 10^-6;
        max_m = [];
        max_avg_turn = 5;
        normalize_data = 1;
        pause_bool = 0;
        
        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,x,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;
        
        
        % Example 3: 3 spirals
    case 3
        n0 = 1000;
        v = (1:5/(n0-1):6)';
        eps = .2;
        rng(1);
        x1 = [v.*cos(v),v.*sin(v)] + .1*sqrt(3)/4*repmat([1,-1],n0,1) + eps*randn(n0,2);
        M = [cos(2*pi/3),-sin(2*pi/3);sin(2*pi/3),cos(2*pi/3)];
        x2 = (M*x1')' + .1*repmat([0,sqrt(3)/4],n0,1) + eps*randn(n0,2);
        x3 = ((M*M)*x1')' + .1*sqrt(3)/4*repmat([-1,-1],n0,1) + eps*randn(n0,2);
        x = [x1;x2;x3];
        n = 3*n0;
        mass = ones(1,n)./n;
        
        x_mean = mass*x;
        x_var = mass*sum((x-repmat(x_mean,n,1)).^2,2);
        x = (x-repmat(x_mean,n,1))/sqrt(x_var); % normalized x
        coeff = pca(x);
        m = 20;
        y0 = bsxfun(@times,(-1.5:3.25/(m-1):1.75)',repmat(coeff(:,1)',m,1));
        cut_indices0 = [];
        
        lambda1 = .002;
        lambda2 = .33;
        
        rho = .05;
        tol = 10^-3;
        max_m = 100;
        max_avg_turn = 10;
        normalize_data = 0;
        pause_bool = 0;
        
        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,x,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;
        
        
        % Example 4: Arc with spiral around it
    case 4
        n1 = 3000;
        n2 = 1250;
        v1 = 0:0.8*pi/n1:0.8*pi';
        v2 = 0:0.8*pi/n2:0.8*pi';
        n1=n1+1; n2=n2+1;
        
        rng(1);
        h1 = .06;
        h2 = .06;
        r1 = h1*unifrnd(-1,1,n1,3);
        r2 = h2*unifrnd(-1,1,n2,3);
        x1 = [ (1+0.2*cos(12*v1)).*cos(v1);(1+0.2*cos(12*v1)).*sin(v1);0.2*sin(12*v1)  ]'+ r1;
        x2 = [ cos(v2);sin(v2);0*v2  ]'+r2;
        x = [x1;x2];
        n = n1+n2;
        mass = 1/n*ones(1,n);
        
        %k = 30; y0 = datasample(x,k); cut_indices0 = (1:k-1)';
        y0 = []; cut_indices0 = [];
        
        crit_dens = .075;
        lambda1 = .0006;
        lambda2 = 4/3*sqrt(lambda1/crit_dens);
        
        rho = 1;
        tol = 10^-4;
        max_m = [];
        max_avg_turn = 15;
        normalize_data = 0;
        pause_bool = 0;
        
        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,x,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;
        
        
        % Two arcs
    case 5
        n1 = 1000;
        v = (0:pi/(n1-1):pi)';
        x1 = [cos(v),sin(v)];
        r = .8;
        n2 = floor(r*n1);
        v = (0:pi/(n2-1):pi)';
        x2 = r*[cos(v),sin(v)];
        n = n1+n2;
        h = .04;
        rng(2);
        x = [x1;x2] + h*randn(n,2);
        mass = 1/n*ones(1,n);
        
        m = 30;
        y0 = [(-1:2/(m-1):1)',zeros(m,1)];
        cut_indices0 = [];
        
        crit_dens = .8/(2*pi);
        lambda1 = h^2*4/(1.8*pi);
        lambda2 = 4/3*sqrt(lambda1/crit_dens);
        
        rho = [];
        tol = 10^-5;
        max_m = [];
        max_avg_turn = 15;
        normalize_data = 0;
        pause_bool = 0;
        
        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,x,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;
        
        
        % Example 6: A noisy NxN grid with background noise
    case 6
        d = 3;
        N = 2;
        L = N+1;
        n = L*200; %points per line
        eps = .05;
        n1 = 2*n*N; %number of points from signal + noise
        n2 = n1; %number of background clutter points
        x = zeros(n1+n2,d);
        mass = 1/(n1+n2)*ones(1,n1+n2);
        rng(1);
        
        for i=1:N
            % vertical data
            x1(:,1) = i*ones(n,1) + eps*randn(n,1);
            x1(:,2) = L/n:L/n:L;
            x1(:,3:d) = eps*randn(n,d-2);
            % horizontal data
            x2(:,1) = L/n:L/n:L;
            x2(:,2) = i*ones(n,1) + eps*randn(n,1);
            x2(:,3:d) = eps*randn(n,d-2);
            x(2*(i-1)*n+1:2*i*n,:) = [x1;x2];
        end
        
        x3 = [rand(n2,2)*L,zeros(n2,d-2)];
        x3 = [rand(n2,2)*L,rand(n2,d-2)*L/2-L/4];
        x(n1+1:n1+n2,:) = x3;
        
        %k = 100; y0 = datasample(x,k); cut_indices0 = (1:k-1)';
        y0 = []; cut_indices0 = [];
        
        lambda1 = .0007;
        lambda2 = .2;
        
        rho = .05;
        tol = 10^-4;
        max_m = 50;
        max_avg_turn = [];
        normalize_data = 1;
        pause_bool = 0;
        
        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,x,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;
        
        
        %Example 7: Cosmic data
    case 7
        x = load('data/cosmic_good.mat');
        x = x.cosmic_good(:,1:2);
        [n,~] = size(x);
        mass = 1/n*ones(1,n);
        
        k = 40; y0 = datasample(x,k); cut_indices0 = (1:k-1)';
        
        lambda1 = .00075;
        lambda2 = .18;
        
        rho = .9;
        tol = 10^-5;
        max_m = [];
        max_avg_turn = 15;
        normalize_data = 1;
        pause_bool = 0;
        
        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,x,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;
        
        %Example 8: Filamentary 3-d data
    case 8
        data = load('data/points3d.mat');
        x = data.v;
        n = length(x(:,1));
        mass = 1/n*ones(1,n);
        
        rng(1);
        y0 = []; cut_indices0 = [];
        %k = 30; y0 = datasample(x,k); cut_indices0 = (1:k-1)';
        
        lambda1 = .00004;
        lambda2 = .08;
        
        rho = .01;
        tol = 10^-4;
        max_m = [];
        max_avg_turn = 40;
        normalize_data = 1;
        pause_bool = 0;
        
        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,x,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;
        
        %Example 9: 'MPPC' letters with background noise
    case 9
        rng(1);
        n1 = 3000;
        
        v_2 = 0:1/(.1*n1-1):1;
        x1_1 = [zeros(.1*n1,1),v_2']; % first segment of 'M'
        v1 = 0:1/2/(.05*n1-1):1/2;
        v_2 = 1:-(1/3)/(.05*n1-1):2/3;
        x1_2 = [v1',v_2']; % second segment of 'M'
        v1 = 1/2:1/2/(.05*n1-1):1;
        v_2 = 2/3:(1/3)/(.05*n1-1):1;
        x1_3 = [v1',v_2']; % third segment of 'M'
        v_2 = 1:-1/(.1*n1-1):0;
        x1_4 = [ones(.1*n1,1),v_2']; % fourth segment of 'M'
        x1 = [x1_1;x1_2;x1_3;x1_4]; %'M'
        
        x1_1 = [1.75*ones(.1*n1,1),(0:1/(.1*n1-1):1)'];
        v1 = -3*pi/8:7/8*pi/(.125*n1-1):pi/2;
        x1_2 = [1.75+.5*cos(v1)',.7+.3*sin(v1)']; % first 'P'
        x2 = [x1_1;x1_2];
        
        x1_1 = [3*ones(.1*n1,1),(0:1/(.1*n1-1):1)'];
        v1 = -3*pi/8:7/8*pi/(.125*n1-1):pi/2;
        x1_2 = [3+.5*cos(v1)',.7+.3*sin(v1)']; % second 'P'
        x3 = [x1_1;x1_2];
        
        v1 = pi/4:3*pi/2/(.25*n1-1):7*pi/4;
        x4 = [4.75+.5*cos(v1)',.5+.5*sin(v1)']; % 'C"
        x = [x1;x2;x3;x4] + .03*randn(n1,2);
        
        n0 = 6000; % number of background clutter
        x0 = [unifrnd(-.75,5.75,n0,1),unifrnd(-.5,1.5,n0,1)];
        x = [x;x0];
        n = n0 + n1;
        mass = ones(1,n)/n;
        
        %k = 50; y0 = datasample(x,k); cut_indices0 = (1:k-1)';
        y0 = []; cut_indices0 = [];
        
        lambda1 = .0007;
        crit_dens = .029;
        lambda2 = 4/3*sqrt(lambda1/crit_dens);
        
        rho = .01;
        tol = 10^-4;
        max_m = [];
        max_avg_turn = 10;
        normalize_data = 0;
        pause_bool = 0;
        
        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,x,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;
        
        
        %Example 10: Zebrafish images
    case 10
        data = load('data/zfishf_51_170.mat'); %289x321x3
        x = data.x;
        clear data;
        n = length(x(:,1));
        mass = 1/n*ones(1,n);
        y0 = []; cut_indices0 = [];
        
        lambda1 = .001;
        lambda2 = 5;
        
        rho = 1;
        tol = 10^-2;
        max_m = [];
        max_avg_turn = 30;
        normalize_data = 1;
        pause_bool = 0;
        
        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,x,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;
        m = length(yfinal(:,1));
        Y = zeros(100,100,m);
        for i=1:m
            Y(:,:,i) = reshape(yfinal(i,:),100,100);
        end
        Y = uint8(Y);
        implay(Y);


        
end

end

