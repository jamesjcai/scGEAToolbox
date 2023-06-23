
C1=quantumCircuit(ryGate(1,0.5*pi))
C2=C1;
isequal(simulate(C1),simulate(C2))

C2.Gates.Angles=0.25*pi;
isequal(simulate(C1),simulate(C2))

% https://www.mathworks.com/matlabcentral/answers/1987423-change-of-parameters-of-gates-in-quantumcircuit-does-not-take-effect