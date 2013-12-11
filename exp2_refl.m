% This script reproduces figure 5 from the paper
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

% inversion
fh  = @(x)misfit(x,Q,D,model,true,false);
fhc = @(x)misfit(x,Q,D,model,true,true);

options.maxiter = 1000;
options.tol     = 1e-3;
options.M       = 5;
options.fid     = 1;
options.cgtol   = 1e-1;
options.cgmax   = 200;
options.write   = 1;

tic
options.method  = 'LBFGS';
[m1,info1]      = GQN(fh,m0,options);
T1    = toc;
cost1 = 2*info1(:,2)+ 3*cumsum(info1(:,6));
for k = 1:size(info1,1)
    mk = dlmread(['x_' num2str(k-1) '.dat']);
    error1(k) = norm(m - mk)/norm(m - m0);
end
delete('x*.dat');
display(['Elapsed time: ' num2str(T1) ' s,  PDE solves: ' num2str(cost1(end))]);


tic
options.method  = 'GN';
[m2,info2]      = GQN(fh,m0,options);
T2 = toc;
cost2 = 2*info2(:,2)+ 3*cumsum(info2(:,6));
for k = 1:size(info2,1)
    mk = dlmread(['x_' num2str(k-1) '.dat']);
    error2(k) = norm(m - mk)/norm(m - m0);
end
delete('x*.dat');
display(['Elapsed time: ' num2str(T2) ' s,  PDE solves: ' num2str(cost2(end))]);


tic
options.method  = 'GN';
[m3,info3]      = GQN(fhc,m0,options);
T3 = toc;
cost3 = 2*info3(:,2)+ 3*cumsum(info3(:,6));
for k = 1:size(info3,1)
    mk = dlmread(['x_' num2str(k-1) '.dat']);
    error3(k) = norm(m - mk)/norm(m - m0);
end
delete('x*.dat');
display(['Elapsed time: ' num2str(T3) ' s,  PDE solves: ' num2str(cost3(end))]);


% plot
figure;
imagesc(1e-3*x,1e-3*z,reshape(1./sqrt(m),n)); axis equal tight;caxis([min(1./sqrt(m)) max(1./sqrt(m))]);
title('ground truth [km/s]','fontsize',20);colorbar;
set(gca,'fontsize',20);xlabel('x [km]','fontsize',20);ylabel('z [km]','fontsize',20);

figure;
imagesc(1e-3*x,1e-3*z,reshape(1./sqrt(m1),n));axis equal tight;caxis([min(1./sqrt(m)) max(1./sqrt(m))]); title('LBFGS','fontsize',20);
set(gca,'fontsize',20);xlabel('x [km]','fontsize',20);ylabel('z [km]','fontsize',20);

figure;
imagesc(1e-3*x,1e-3*z,reshape(1./sqrt(m2),n)); axis equal tight;caxis([min(1./sqrt(m)) max(1./sqrt(m))]);title('GN w/o correction','fontsize',20);
set(gca,'fontsize',20);xlabel('x [km]','fontsize',20);ylabel('z [km]','fontsize',20);

figure;
imagesc(1e-3*x,1e-3*z,reshape(1./sqrt(m3),n)); axis equal tight;caxis([min(1./sqrt(m)) max(1./sqrt(m))]);title('GN w correction','fontsize',20);
set(gca,'fontsize',20);xlabel('x [km]','fontsize',20);ylabel('z [km]','fontsize',20);

figure;
loglog(1+info1(:,1),info1(:,4),'k',1+info2(:,1),info2(:,4),'r',1+info3(:,1),info3(:,4),'b','linewidth',2);
legend('LBFGS','GN w/o corr','GN w corr');
set(gca,'fontsize',20);xlabel('iteration','fontsize',20);ylabel('misfit','fontsize',20);

figure;
loglog(1+info1(:,1),info1(:,5),'k',1+info2(:,1),info2(:,5),'r',1+info3(:,1),info3(:,5),'b','linewidth',2);
legend('LBFGS','GN w/o corr','GN w corr');
set(gca,'fontsize',20);xlabel('iteration','fontsize',20);ylabel('||g||_2','fontsize',20);

figure;
semilogx(1+info1(:,1),error1,'k',1+info2(:,1),error2,'r',1+info3(:,1),error3,'b','linewidth',2);
legend('LBFGS','GN w/o corr','GN w corr');
set(gca,'fontsize',20);xlabel('iteration','fontsize',20);ylabel('reconstruction error','fontsize',20);

for k = 1:7
    print(k,'-depsc',[mfilename '_' num2str(k)]);
end

