% to get the data for this plot:
%   hamilton_harm.py 10 10.0 
%	produces the potential data and the input to the next program
%   power_ham.py 10 0 0 invertHarmonic_L10_xmax10.0.json evecHarm_L10_xmax10.0
%       produces the eigenvector data (currently with maxnorm 1, thus the
%       scale factor below to normalize the integral to 1)
%
x = linspace(-10,10,1024);
y = 0.5*ones(1,1024);
potential = load('harmonicPot_L10_xmax10.0')';
evec0 = load('evecHarm_L10_xmax10.0')';
% normalize eigenvector so that integral |evec0|^2 = 1
delta_x = 20/1023.0;
norm = sum(evec0.*evec0)*delta_x;
evec0 = evec0/sqrt(norm);
% make plot
plot(x,potential)
hold
plot(x,y,'r--')
plot(x,5*evec0+y,'red')
axis([-5,5,0,20])
hold
% for comparison purposes, here is the analytic form for the wavefunction
analytic = exp(-x.*x/2);
norm1 = sum(analytic.*analytic)*delta_x;
analytic = analytic/sqrt(norm1);
