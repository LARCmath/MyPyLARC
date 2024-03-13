%\*****************************************************************
%* Copyright (C) 2014-2024, Institute for Defense Analyses        *
%* 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           *
%* This material may be reproduced by or for the US Government    *
%* pursuant to the copyright license under the clauses at DFARS   *
%* 252.227-7013 and 252.227-7014.                                 *
%*                                                                *
%* LARC (Linear Algebra via Recursive Compression)                *
%* Authors:                                                       *
%*   - Steve Cuccaro (IDA-CCS)                                    *
%*   - John Daly (LPS)                                            *
%*   - John Gilbert (UCSB, IDA adjunct)                           *
%*   - Mark Pleszkoch (IDA-CCS)                                   *
%*   - Jenny Zito (IDA-CCS)                                       *
%*                                                                *
%* Additional contributors are listed in "LARCcontributors".      *
%*                                                                *
%* Questions: larc@super.org                                      *
%*                                                                *
%* All rights reserved.                                           *
%*                                                                *
%* Redistribution and use in source and binary forms, with or     *
%* without modification, are permitted provided that the          *
%* following conditions are met:                                  *
%*   - Redistribution of source code must retain the above        *
%*     copyright notice, this list of conditions and the          *
%*     following disclaimer.                                      *
%*   - Redistribution in binary form must reproduce the above     *
%*     copyright notice, this list of conditions and the          *
%*     following disclaimer in the documentation and/or other     *
%*     materials provided with the distribution.                  *
%*   - Neither the name of the copyright holder nor the names of  *
%*     its contributors may be used to endorse or promote         *
%*     products derived from this software without specific prior *
%*     written permission.                                        *
%*                                                                *
%* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         *
%* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
%* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       *
%* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       *
%* DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        *
%* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   *
%* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   *
%* NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   *
%* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       *
%* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      *
%* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   *
%* OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, *
%* EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             *
%*                                                                *
%*****************************************************************/

% to get the data for this plot:
%   hamilton_morse.py 10 5.0 
%	produces the potential data and the input to the next program
%   power_ham.py 10 0 0 invertMorse_L10_rmax5.0.json evecMorse_L10_rmax5.0
%       produces the eigenvector data (currently with maxnorm 1, thus the
%       scale factor below to normalize the integral to 1)
%

% Set up Morse potential parameters, calculate eigenenergy of ground state
re = 1.28;      	% distance for miminum energy in potential (Å)
de = 4.17;		% well depth (eV)
alpha = 1.85;   	% parameter determining well width (Å^{-1})
hbar = 6.582119e-16;	% Planck's constant/2\pi (eV⋅s)
c = 2.99792458e18;	% speed of light (Å/s)
m = 9.15203813e8;	% reduced mass of system (eV/c^2)
% * 1+\floor(lambda-0.5) is # of bound states in system (unitless)
lambda = sqrt(2*m*de)/alpha/hbar/c;
% * E_0 is the ground state energy (eV)
C_0 = power(hbar*alpha,2)/2/m*c*c;
E_0 = -C_0*power(lambda-0.5,2);
% need scaled variables for analytic eigenvector calculation
r = linspace(0,5,1024);
x = alpha*r;
xe = alpha*re;
% load calculated values
potential = load('morsePot_L10_rmax5.0')';
evec0 = load('evecMorse_L10_rmax5.0')';
% normalize so that integral |evec0|^2 = 1
delta_r = 5.0/1023;
norm = sum(evec0.*evec0)*delta_r;
evec0 = evec0/sqrt(norm);
% calculate analytic form values
% * z = 2\lambda\exp{-x+x_e}
z = 2*lambda*exp(-x+xe);
% * N_0 = \left[\frac{0!(2\lambda-1)}{\Gamma(2\lambda)}\right]^{1/2}
N0 = sqrt((2*lambda-1)/gamma(2*lambda));
% * \Psi_0(z) = N_0 z^{\lambda-0.5}e^{-z/2}
%	(the term L_0^{y}(z) == 1 is left out)
analytic = N0*power(z,lambda-0.5).*exp(-z/2);
% normalize analytic so that integral |analytic|^2 = 1
norm1 = sum(analytic.*analytic)*delta_r;
analytic = analytic/sqrt(norm1);
% make plot
eval0 = -3.997578381116;
y = eval0*ones(1,1024);
plot(r,potential,'k')
hold
plot(r,y,'m--')
plot(r,evec0+y,'m')
axis([0,5,-5,2])
% format the plot labels
xstr = sprintf('distance between particles (%c)',197);
xlabel(xstr)
ylabel('potential (eV) / state amplitude')
titlestr1 = sprintf('Ground state of Morse potential for ^{40}Ar-p interaction');
titlestr2 = sprintf('D_e=4.17eV, r_e=1.28%c, α=1.85%c^{-1}',197,197);
title({titlestr1, titlestr2})
% make fonts big enough to read
set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(findall(gcf,'-property','FontName'),'FontName','Helvetica')
hold
% When you change the font size, this doesn't work right - you have to manually
% stretch the figure window until all the text is visible, then do a manual
% save, then re-stretch till the figure looks right, then save again.
%print('morseplot.jpg','-djpeg')

