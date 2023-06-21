% cd 'D:\GitHub\scGEAToolbox\+run\external\py_SERGIO'

% Step 1: Draw GRN
% gui.graph_gui

% Step 2: Save GRN as A
load("example_grn.mat");

% Setp 3: Write A to regs.txt and targets.txt
pkg_e_writesergiogrn(A)

% Step 4: Run SERGIO
X=run.py_SERGIO(A,5000);


%%
B=sc_pcnet(X);

figure;
subplot(2,2,1)
p=plot(digraph(A)); title('Ground Truth')

subplot(2,2,2)
i_plotit(B,0.6,p)

subplot(2,2,3)
i_plotit(B,0.7,p)

subplot(2,2,4)
i_plotit(B,0.8,p)


function i_plotit(B,c,p)
    p2=plot(digraph(ten.e_filtadjc(B,c))); 
    p2.XData=p.XData; 
    p2.YData=p.YData; 
    title(sprintf('PCNet cutoff=%g',c));
end
