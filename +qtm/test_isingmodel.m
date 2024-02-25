clear all
close all 
Nspin=49; % number of spins. Must be square of integer 
J=2; %  constant of interaction 
h=0; % external magnetic field 
Esi=10; % energy of system (microcanonic distribution) 
NTrial=1000; % number of trials 
[Es,Ed,SpM,A,S] = qtm.Ising2(Nspin,J,h,Esi,NTrial); 
% plot of energy of system 
% plot of energy of daemon 
[~,m]=size(Ed); % [~,m], that means that you just want the second output of your function, and do not care the first one. 
% is equal to Nspin*Ntrial and considered as time 
i=1:m; 

figure(1); plot(i,Es(i)) 
figure(2); plot(i,Ed(i)) 

% plot of spin coniguration 
m1=floor(m/3); 
m2=2*m1; 
m3=3*m1; 
i=1:sqrt(Nspin); 
j=1:sqrt(Nspin); 
V(i,j)=0; 
U(i,j)=0; 
Z(i,j)=0; 

figure 
subplot(4,1,1);quiver3(Z,U,V,S(:,:,1)); colormap white; 
title('Spins initial configuration') 
subplot(4,1,2);quiver3(Z,U,V,S(:,:,m1)); colormap white 
title(['Spins configuration for time ', num2str(m1)]) 
subplot(4,1,3); quiver3(Z,U,V,S(:,:,m2)); colormap white 
title(['Spins configuration for time ', num2str(m2)]) 
subplot(1,1,1); quiver3(Z,U,V,S(:,:,m)); colormap white 
 title(['Spins configuration for time ', num2str(m)]) 

mean(Es) % mean energy of system
