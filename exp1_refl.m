% This script reproduces figure 2(b) from the paper
%
% -------------------------------------------------------------------------
%      Copyright (C) 2013 Tristan van Leeuwen
%                         Centrum Wiskunde & Informatica
%                         Tristan.van.Leeuwen@cwi.nl
%  
%      This program is free software: you can redistribute it and/or modify
%      it under the terms of the GNU General Public License as published by
%      the Free Software Foundation, either version 3 of the License, or
%      (at your option) any later version.
%  
%      This program is distributed in the hope that it will be useful,
%      but WITHOUT ANY WARRANTY; without even the implied warranty of
%      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%      GNU General Public License for more details.
%  
%      You should have received a copy of the GNU General Public License
%      along with this program.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

% define velocity model
n  = [51 51];
h  = [20 20];
z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]*h(2);

v0 = 2000*ones(n);
dv = zeros(n);
dv(11:21,11:21) = 100;
dv(31:41,31:41) = -100;
m  = 1e6./(v0(:) + dv(:)).^2;
m0 = 1e6./(v0(:)).^2;

% grid
model.h = h;
model.n = n;

% frequency
model.f  =  [2.5 5 10];

% acquisition
model.Is = {2,2:50};
model.Ir = {2,2:50};

% sources
Q = speye(length(model.Is{1})*length(model.Is{2}));

% data
D = F(m,Q,model);

% multiply with random amplitude factor
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
w = 5*abs(randn(1,size(D,2),size(D,3)));
D = bsxfun(@times,w,D);

% misfit
fh  = @(x)misfit(x,Q,D,model,true,false);
fhc = @(x)misfit(x,Q,D,model,true,true);

% perturbation
[zz,xx] = ndgrid(z,x);
mr      = m0;
alpha   = linspace(-.2,.2,20);


[fr,gr,Hr] = fh(mr);
%dm         = exp(-1e-4*((xx(:) - mean(x)).^2 + (zz(:) - mean(z)).^2))*mean(m0);
dm         = gr/norm(gr);
%dm         = m - m0;

hr         = Hr(dm);

[frc,grc,Hrc] = fhc(mr);
hrc           = Hrc(dm);

for k = 1:length(alpha)
    f(k)  = fh(mr + alpha(k)*dm);
    q1(k) = fr + alpha(k)*dm'*gr  + .5*alpha(k)^2*dm'*hr;
    q2(k) = fr + alpha(k)*dm'*grc + .5*alpha(k)^2*dm'*hrc;
end

% plot
figure;
plot(alpha,f,'k--',alpha,q1,'k',alpha,q2,'r','linewidth',2); 
legend('true','w/o correction','w correction','location','NorthWest');
xlabel('\alpha','fontsize',20);ylabel('misfit','fontsize',20);
set(gca,'fontsize',20);

% save
print(1,'-depsc',mfilename);
