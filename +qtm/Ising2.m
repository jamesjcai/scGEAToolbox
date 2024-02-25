function [Es,Ed,SpM,A,S] = Ising2(Nspin,J,h,Esi,NTrial)
% input parameters are described in square_grid1.m
% output parameters
% Es is instantaneous energy of system after each step of walk in each trial
%SpM,A,S are spin configuration
Ns=Nspin.^0.5; % number of spin for one side of square
s=ones(Ns,Ns); % initial spin configuration. All spin projections are 1
Esystem=-(J+h)*Nspin; % initial system energy
Edemon=4*J*floor((Esi-Esystem)/(4*J)); % initial energy of daemon.
% on each trial daemon energy is compared with energy of all spin knots of grid
% and energy of knot and daemon is changed not changed according some
% condition
Es(1)=Esystem; 
Ed(1)=Edemon;
S=s;
k=1;
for i=1:NTrial
  Accept=0;
  for j=1:Nspin
     % random choice of knot
     Ix=floor(Ns*rand(1)+1); 
     Iy=floor(Ns*rand(1)+1);
     % border conditions
     if Ix==1
        Left=Ns;
     else
        Left=Ix-1;
     end

     if Ix==Ns
       Right=1;
     else
        Right=Ix+1;
     end

     if Iy==1
        Down=Ns;
     else  
        Down=Iy-1;
     end

     if Iy==Ns
        Up=1;
     else
        Up=Iy+1;
     end

     % energy change
     de=2*s(Iy,Ix)*(-h+J*(s(Iy,Left)+s(Iy,Right)+s(Down,Ix)+s(Up,Ix)));
     if de<=Edemon % energy change accepted
        s(Iy,Ix)=-s(Iy,Ix);
        Accept=Accept+1;
        Edemon=Edemon-de;
        Esystem=Esystem+de;
     end

     k=k+1;
     Es(k)=Esystem;
     Ed(k)=Edemon;
     A(k-1)=Accept;
     s1=sum(s);
     SpM(k)=sum(s1);
     S=cat(3,S,s);
  end
end
A=A/NTrial;
end