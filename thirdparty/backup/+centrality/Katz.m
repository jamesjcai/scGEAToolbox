function x=Katz(adj)
x = (eye(n)-(beta/max_eig)*adj)\ones(n,1);

% http://guettel.com/rktoolbox/examples/html/example_networks.html
