function [ y_new,cut_indices,I,iter ] = mppc( y0,cut_indices,x,mass,lambda1,lambda2,tol,rho0,...
    max_m,max_avg_turn,normalize_data,pause_bool,plotting)
%MPPC the main function for the algorithm that finds local minimizers of
%the MPPC functional, as introduced in 
% Kirov, S. & Slep?ev, D., Multiple Penalized Principal Curves: Analysis 
% and Computation, J Math Imaging Vis (2017). doi:10.1007/s10851-017-0730-8

% INPUTS:
% y0 - the initial polygonal curve guess. Dimension is mxd, where m is the
% number of points on the curve and d is the space dimension.
%
% cut_indices - empty vector if y0 is connected. If not, this column vector
% contains the indices in y0 after which a disconnection occurs.
%
% x - the data. Dimension is nxd, where n is the number of data points, and
% d is the space dimension. The data can be normalized (see below).
%
% mass - a 1xn row vector signifying the mass of each data point. Usually
% each entry = 1/n.
%
% lambda1 - the coefficient for length penalty in the objective functional.
%
% lambda2 - the coefficient for the penalty on the number of curves.
%
% The following are all optional inputs: pass empty vectors for default
% values.
%
% tol - stopping criteria for convergence, based on relative difference of
% successive energies. The algorithm terminates when 
% old energy - new energy < tol*(new energy). Default value is 10^-3.
%
% rho0 - initial parameter for the split bregman/ADMM algorithm. Lower values
% correspond to larger "step sizes" of y. rho may change throughout the
% algorithm based on the primal and dual residuals. Default value = lambda1
%
% max_m - maximum desired number of points on curve. It may end up slightly
% higher due to no singletons left behind when cutting consecutive edges.
%
% max_avg_turn - points will continue to be added to y until m>max_m or the
% average turn degree < max_avg_turn.
%
% normalize_data - will standard-normalize data iff set to 1/true.
%
% pause_bool - will pause plots after each iteration iff set to 1/true.
%
% OUTPUTS:
% y_new - the found curve(s) stored as a mxd matrix, where m is the number
% of points on the curve(s), and d is the dimension.
%
% cut_indices - a vector containing the indices of y_new after which a new
% curve begins.
%
% I - a list of 1xn assignments containg which of the 1,...,m indices each
% of the n data points is closest to.
%
% differs from public version in that here modify spacing has no effect if
% it does not decrease the energy. In the public version, it is given a
% chance. 

MAX_ITER = 900; %maximum number of iterations
MAX_INNER_ITER = 500; %maximum number of inner ADMM iterations
r = 3; %param for modifySpacing
cut_freq = 10; %frequency with which to check topology
respace_freq = 5; %frequency with which to check parametrization
% plotting = 0; %plot if and only if = 1

if isempty(tol), tol = 10^-3; end
num_d = ceil(-log(tol)/log(10)); %number of decimals to display of energy
if isempty(rho0)
    rho = lambda1;
else
    rho = rho0;
end

assert(isrow(mass),'mass should be a 1xn row vector');
assert(isempty(cut_indices) | iscolumn(cut_indices), 'cut_indices should be a column vector (or empty)');
if normalize_data 
    n = length(x(:,1));
    x_mean = mass*x;
    x_var = mass*sum((x-repmat(x_mean,n,1)).^2,2);
    x = (x-repmat(x_mean,n,1))/sqrt(x_var);
end
if isempty(y0)
    fprintf('Initializing...'); tic;
    [y0,cut_indices,~] = initialize(x,mass,lambda1,lambda2);
    toc;
else 
    if normalize_data
        m = length(y0(:,1));
        y0 = (y0 - repmat(x_mean,m,1))/sqrt(x_var);
    end
end
y = y0;
if plotting, plotADP(x,y0,cut_indices,[],[],0,[0,.8,.8],'-'); end

eo = 0; % initial boolean for respacing of even/odd points
iter = 0;

y_new = y;
z_new = y_new(2:end,:)-y_new(1:end-1,:) ; b = zeros(size(z_new));
I = dsearchn(y_new,x)';
[energy,cont_energy] = calculateEnergy(y_new,x,mass,lambda1,lambda2,I,cut_indices);
num_comps = length(cut_indices)+1;

m = length(y_new(:,1));
rgb_c = [0,.8,0]; ls = '-'; %color and linestyle of curve when plotting
if pause_bool, plotADP(x,y_new,cut_indices,[],[],0,rgb_c,ls); pause; end;

fprintf(['\n iter = %d       E = %8.',num2str(num_d),'E      cont_E = %8.',...
    num2str(num_d),'E      k = %d     m = %d'],iter,energy,cont_energy,num_comps,m);

check_cuts = 1; check_conn = 1; check_space = 1; check_singletons = 1;
more_prec = 0; energy_prev = inf; inner_iter = 1; 
avg_turn = [];

while (~isequal(size(y),size(y_new)) || energy_prev-energy>tol*energy || energy_prev<energy ...
        || check_cuts || check_conn || check_space || check_singletons || inner_iter==1) && iter<MAX_ITER
    
    iter = iter+1; y = y_new; z = z_new; energy_prev = energy;
    if pause_bool, fprintf('\n press a key to bregman step'), pause; end
    
    [y_new,z_new,b,I,energy,cut_indices] = updateyzb(y_new,z_new,b,x,mass,lambda1,rho,cut_indices,I);
    
    if length(z_new(:,1))==length(z(:,1)) && more_prec && ~isempty(z_new) %Check residuals to determine whether to change rho
        [r_pri,r_dual] = computeResiduals(y_new,z_new,z,rho);
        if r_pri>5*r_dual, rho = rho*1.3;
            if ~isempty(rho0), fprintf('    rho = %3.2E',rho); end
        end
        if 5*r_pri<r_dual, rho = rho/2;
            if ~isempty(rho0), fprintf('    rho = %3.2E',rho); end
        end
    end

    num_comps = length(cut_indices)+1;
    energy = energy + lambda1*lambda2*(num_comps-1);
    energy0 = 0;
    inner_iter = 1;
    while (energy>energy_prev || (more_prec && (energy0-energy>tol*energy || energy0-energy<0))) ...
            && abs(energy-energy0)>10^-14 && inner_iter<MAX_INNER_ITER %do inner ADMM iterations
        energy0 = energy;
        z = z_new;
        [y_new,z_new,b,I,energy,cut_indices] = updateyzb(y_new,z_new,b,x,mass,lambda1,rho,cut_indices,I);
        energy = energy + lambda1*lambda2*(length(cut_indices));
        inner_iter = inner_iter+1;
        if pause_bool, fprintf(['\n inner iter = %d       energy = %8.',num2str(num_d),'E '],inner_iter,energy);
            plotADP(x,y_new,cut_indices,mass,I,0,rgb_c,ls), pause;
        end
    end
    
    if inner_iter > 1 && ~isempty(z_new) %Check residuals to determine whether to change rho
        [r_pri,r_dual] = computeResiduals(y_new,z_new,z,rho);
        if r_pri>5*r_dual, rho = rho*1.3;
            if ~isempty(rho0), fprintf('    rho = %3.2E',rho); end
        end
        if 5*r_pri<r_dual, rho = rho/2;
            if ~isempty(rho0), fprintf('    rho = %3.2E',rho); end
        end
        fprintf('\n %d inner ADMM iterations',inner_iter);
        if more_prec == 1
            fprintf('  (more precision in minimization obtained)');
            more_prec = 0;
        end
    end
    
    [m,~] = size(y_new);
    I = dsearchn(y_new,x)';
    [reproj_energy, cont_energy] = calculateEnergy(y_new,x,mass,lambda1,lambda2,I,cut_indices);
    if energy-reproj_energy < 10*tol*energy
        more_prec = 1;
    end
    energy = reproj_energy;
    fprintf(['\n iter = %d       E = %8.',num2str(num_d),'E      cont_E = %8.',...
        num2str(num_d),'E      k = %d     m = %d'],iter,energy,cont_energy,num_comps,m);
    if ~isempty(avg_turn)
        fprintf('    avg turn = %.2f',avg_turn);
    end
    
    if plotting, plotADP(x,y_new,cut_indices,mass,I,0,rgb_c,ls); drawnow; end
    
    if mod(iter+4,cut_freq) == 0 %check cutting
        if pause_bool, fprintf('\n press a key to cut edges'), pause; end
        check_cuts = 1;
        [energy, ~] = calculateEnergy(y_new,x,mass,lambda1,lambda2,I,cut_indices);
        [y_new,z_new,b,I,cut_indices_new] = energyCut(x,mass,y_new,z_new,b,[],cut_indices,lambda1,lambda2);
        if isequal(cut_indices,cut_indices_new)
            check_cuts = 0;
        else
            cut_indices = cut_indices_new; more_prec = 0;
            [energy, cont_energy] = calculateEnergy(y_new,x,mass,lambda1,lambda2,I,cut_indices);
            fprintf(['\n cutting:    %s   E = %8.',num2str(num_d),'E      cont_E = %8.',...
                num2str(num_d),'E      k = %d'],repmat(' ',1,floor(log(iter)/log(10))),energy,cont_energy,length(cut_indices)+1);
        end
        if pause_bool, plotADP(x,y_new,cut_indices,mass,I,0,rgb_c,ls); pause; end;
    end
    
    if mod(iter+3,cut_freq)==0 % check singletons
        if pause_bool, fprintf('\n press a key to check singletons'), pause; end
        [y_new,z_new,b,cut_indices_new] = checkSingletons(x,mass,I,y_new,z_new,b,cut_indices,lambda1,lambda2);
        if ~isequal(cut_indices,cut_indices_new)
            check_singletons = 1;
            cut_indices = cut_indices_new;
            I = dsearchn(y_new,x)'; more_prec = 0;
            [energy,cont_energy] = calculateEnergy(y_new,x,mass,lambda1,lambda2,I,cut_indices);
            fprintf(['\n singletons:   %s E = %8.',num2str(num_d),'E      cont_E = %8.',...
                num2str(num_d),'E      k = %d'],repmat(' ',1,floor(log(iter)/log(10))),energy,cont_energy,length(cut_indices)+1);
        else
        check_singletons = 0;
        end
        if pause_bool, plotADP(x,y_new,cut_indices,mass,I,0,rgb_c,ls); pause; end;
    end
    
    if mod(iter+2,cut_freq)==0 % check connecting
        if isempty(cut_indices)
            check_conn = 0;
        else
            check_conn = 1;
            if pause_bool, fprintf('\n press a key to connect components'), pause; end
            [y_new,z_new,b,I,~,cut_indices] = remZeroMassPoints(y_new,z_new,b,mass,I,cut_indices);
            [energy, ~] = calculateEnergy(y_new,x,mass,lambda1,lambda2,I,cut_indices);
            [y_new,cut_indices_new] = connectCompsE(x,y_new,mass,I,cut_indices,lambda1,lambda2);
            if ~isequal(cut_indices,cut_indices_new)
                cut_indices = cut_indices_new; check_space = 1;
                I = dsearchn(y_new,x)'; more_prec = 0; 
                z_new = y_new(2:end,:)-y_new(1:end-1,:) ; b = zeros(size(z_new));
                [energy, cont_energy] = calculateEnergy(y_new,x,mass,lambda1,lambda2,I,cut_indices);
                fprintf(['\n connecting: %s   E = %8.',num2str(num_d),'E      cont_E = %8.',...
                    num2str(num_d),'E      k = %d'],repmat(' ',1,floor(log(iter)/log(10))),energy,cont_energy,length(cut_indices)+1);
            else
                check_conn = 0;
            end
            if pause_bool,plotADP(x,y_new,cut_indices,mass,I,0,rgb_c,ls); pause; end;
        end
    end
    
    if mod(iter+1,respace_freq)==0 % check parametrization
        if pause_bool, fprintf('\n press a key to re-space points'), pause; end
        eo = mod(eo+1,2);
        if m<length(x(:,1))
            [num_add,avg_turn] = numPointsAdd(y_new,cut_indices,max_avg_turn,max_m,lambda2);
            check_space = 1;
            if num_add == 0, check_space = 0; end
        else
            num_add = 0;
        end
        y_old = y_new; I_old = I; cut_indices_old = cut_indices;
        energy_old = energy; z_old = z; b_old = b;
        [y_new,~,~,cut_indices,I] = modifySpacingAllComps(y_new,z_new,b,x,mass,I,r,cut_indices,num_add,eo);
        if num_add>0 && length(y_new(:,1)) == length(y(:,1)) %points were not successully added
            check_space = 0;
        end
        z_new = y_new(2:end,:)-y_new(1:end-1,:) ; b = zeros(size(z_new));
        %I = dsearchn(y_new,x)'; 
        [energy,cont_energy] = calculateEnergy(y_new,x,mass,lambda1,lambda2,I,cut_indices);
        if energy>=energy_old
            y_new = y_old; I = I_old; cut_indices = cut_indices_old;
            energy = energy_old;
            z_new = z_old; b = b_old;
            check_space = 0;
        else
            more_prec = 0; num_comps = length(cut_indices)+1;
            fprintf(['\n respacing: %s    E = %8.',num2str(num_d),'E      cont_E = %8.',...
                num2str(num_d),'E                m = %d'],repmat(' ',1,floor(log(iter)/log(10))+floor(log(num_comps)/log(10))),...
                energy,cont_energy,length(y_new(:,1)));
            if pause_bool, plotADP(x,y_new,cut_indices,mass,I,0,rgb_c,ls); pause; end;
        end
    end
end

if normalize_data %un-normalize data
    x = sqrt(x_var)*x + repmat(x_mean,n,1);
    m = length(y_new(:,1));
    y_new = sqrt(x_var)*y_new + repmat(x_mean,m,1);
end

if plotting, plotADP(x,y_new,cut_indices,mass,I,0,rgb_c,ls); end

I = dsearchn(y_new,x)';
fprintf('\n');
end


function [r_pri,r_dual] = computeResiduals(y_new,z_new,z,rho)
r_pri = norm(y_new(2:end,:)-y_new(1:end-1,:)-z_new);
z_diff = z_new - z;
r_dual = rho*norm([z_diff(1,:);z_diff(1:end-1,:)-z_diff(2:end,:);z_diff(end,:)]);
end
