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
