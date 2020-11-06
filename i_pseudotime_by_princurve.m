function [t,xyz1]=i_pseudotime_by_princurve(s,plotit)

if nargin<2, plotit=true; end

pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/MPPC');
addpath(pth);

        n=size(s,1);
        mass = 1/n*ones(1,n);
        y0 = []; cut_indices0 = [];
        
        crit_dens = .075;
        lambda1 = .0006;
        lambda2 = 4/3*sqrt(lambda1/crit_dens);
        
        rho = 1;
        tol = 10^-4;
        max_m = [];
        max_avg_turn = 15;
        normalize_data = 1;
        pause_bool = 0;
        

        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,s,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;

        t=I;
        xyz1=yfinal;


 if plotit
  hold on
  plot3(xyz1(:,1),xyz1(:,2),xyz1(:,3),'.r','linewidth',2);
 end
end
