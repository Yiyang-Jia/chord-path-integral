%aaa = readmatrix('tranptsHigherAccuracyl=1lz=0.txt');
%aaa = readmatrix('tranptsHigherAccuracyl=0p85lz=0.txt');
%aaa = readmatrix('tranptsHigherAccuracyl=1p5lz=0.txt');
%aaa = readmatrix('tranptsHigherAccuracyl=2lz=0.txt');
aaa = readmatrix('tranptsHigherAccuracyl=6lz=0Sorted.txt');


lnz=6;
lzz=0;

aoriginal =aaa;
aaa = sortrows(aaa);
%dlmwrite('tranptsHigherAccuracyl=6lz=0Sorted.txt',aaa)

betajbegin=1/aaa(1,1);
betajend = 1/aaa(end,1);
kappaend = aaa(end, 2);
datasplitpoint =20; %only low temperature point obey linear relation
aaalinear =aaa(1:datasplitpoint,:);
aaaremain = aaa(datasplitpoint+1:end,:);

   [p,S]=polyfit(log(aaalinear(:,1)),log(aaalinear(:,2)),1);
%inter= @(x) interp1(aaa(:,1),aaa(:,2),x,'spline','extrap');

% template = hgload('phasediagramloglog.fig');
% ax=gca;
% ha1 = findobj(template,'Type','Legend');
% 

 pltx2 = linspace(log(1/betajbegin)-2,log(1/betajend),100);
%% lnz=1, lzz=0 
 
 ff = figure;
subplot(1,2,1)
plot(pltx2, p(1)*pltx2+ p(2),log(aaalinear(:,1)),log(aaalinear(:,2)),"o","MarkerSize",4)
hold on
plot(log(aaaremain(:,1)),log(aaaremain(:,2)),"-*","MarkerSize",4)
hold on
plot(log(1/betajend),log(kappaend),"*","Color","black","MarkerSize",4)
xlim([-10 log(1/betajend)/1.5])
ylim([-4.5 log(kappaend)+1])
 set(gca,'FontSize',15)

titlestr1 = 'Phase diagram, log-log scale';
%titlestr1 = append('Phase diagram ($l_{nz}$ =', num2str(lnz,'%.4g'), ', $l_{zz}$ = ', num2str(lzz,'%.4g'), '), log-log scale');

title(titlestr1,'Interpreter','latex','FontSize',20)
% 
% text(log(1/betajend) -7,log(kappaend)-1,'p-spin phase','FontSize', 16)
% text(log(1/betajend)-3,log(kappaend)-2,'SYK phase','FontSize', 16)
text(log(1/betajend) -7,log(kappaend)-1,'quasi-integrable phase','FontSize', 18)
text(log(1/betajend)-3,log(kappaend)-2,'chaotic phase','FontSize', 18)

 xlabel('$\log[1/(\beta \mathbf{J} )]$','Interpreter','latex','FontSize', 20)
 ylabel('$\log \kappa$','Interpreter','latex','FontSize', 20)
 legstr11 = append('$\log \kappa=',num2str(p(1),'%.4f'),' \log(1/\beta \mathbf{J})$ + ', num2str(p(2),'%.4f'));
 legstr12 = 'first-order transition points at low T';
 legstr13 ='other first-order transition points'
 legstr14 = 'second-order transition point'
 legend(legstr11,legstr12,legstr13,legstr14,'Interpreter','latex', 'Location', 'southeast','FontSize',16)
grid on

subplot(1,2,2)
pltx = linspace(0,1/betajend,100);
%plotdat = inter(pltx);
plotdat = exp(p(2))* pltx.^(p(1));
plot(pltx,plotdat, aaalinear(:,1),aaalinear(:,2),"o","MarkerSize",4)
hold on
plot(aaaremain(:,1),aaaremain(:,2),"-*","MarkerSize",4)
hold on
plot(1/betajend,kappaend,"*","Color","black","MarkerSize",4)
set(gca,'FontSize',15)

secondorderPtstring = append( '(', num2str(1/betajend,'%.3g'), ',', num2str(kappaend,'%.3g'), ')');
text(1/betajend-0.005,kappaend-0.025,secondorderPtstring,'FontSize', 12)
% text(0.5/betajend -0.03,0.5*kappaend+0.25,'p-spin phase','FontSize', 16)
% text(0.5/betajend,0.5*kappaend,'SYK phase','FontSize', 16)
text(0.5/betajend -0.025,0.5*kappaend+0.22,'quasi-integrable phase','FontSize', 18)
text(0.5/betajend,0.5*kappaend,'chaotic phase','FontSize', 18)

xlim([0 1.2/betajend])
ylim([0 kappaend*1.5])
titlestr2 = 'Phase diagram, linear scale';
%titlestr2 = append('Phase diagram ($l_{nz}$ =', num2str(lnz,'%.4g'), ', $l_{zz}$ = ', num2str(lzz,'%.4g'),'), linear scale');
title(titlestr2,'Interpreter','latex','FontSize',20);
xlabel('$1/(\beta \mathbf{J})$','Interpreter','latex','FontSize',20);
 ylabel('$\kappa$','Interpreter','latex','FontSize',20);
  legstr21 = append('$\kappa =', num2str( exp(p(2)),'%.4f'),'(1/\beta \mathbf{J})^{', num2str(p(1),'%.4f'),'}$');
 legstr22 = 'first-order transition points at low T';
 legstr23 = 'other first-order transition points';
 legstr24 = 'second-order transition point'
  legend(legstr21,legstr22,legstr23,legstr24,'Interpreter','latex','Location', 'southeast','FontSize',16)

  
grid on

% % xlim([0 1/12])
 %ylim([0 .6])
 %% lnz>1, lzz=0
 ff = figure;
subplot(1,2,1)
plot(pltx2, p(1)*pltx2+ p(2),log(aaalinear(:,1)),log(aaalinear(:,2)),"o","MarkerSize",4)
hold on
plot(log(aaaremain(:,1)),log(aaaremain(:,2)),"-*","MarkerSize",4)
hold on
plot(log(1/betajend),log(kappaend),"*","Color","black","MarkerSize",4)
xlim([-10 log(1/betajend)/1.5])
ylim([-4.5 log(kappaend)+1])
 set(gca,'FontSize',15)

%titlestr1 = 'Phase diagram, log-log scale';
titlestr1 = append('Phase diagram ($l_{nz}$ =', num2str(lnz,'%.4g'), ', $l_{zz}$ = ', num2str(lzz,'%.4g'), '), log-log scale');

title(titlestr1,'Interpreter','latex','FontSize',20)
% 
% text(log(1/betajend) -7,log(kappaend)-1,'p-spin phase','FontSize', 16)
% text(log(1/betajend)-3,log(kappaend)-2,'SYK phase','FontSize', 16)
text(log(1/betajend) -8,log(kappaend)-0.,'quasi-integrable phase','FontSize', 18)
text(log(1/betajend)-3,log(kappaend)-2,'chaotic phase','FontSize', 18)

 xlabel('$\log[1/(\beta \mathbf{J} )]$','Interpreter','latex','FontSize', 20)
 ylabel('$\log \kappa$','Interpreter','latex','FontSize', 20)
 legstr11 = append('$\log \kappa=',num2str(p(1),'%.4f'),' \log(1/\beta \mathbf{J})$ + ', num2str(p(2),'%.4f'));
 legstr12 = 'first-order transition points at low T';
 legstr13 ='other first-order transition points'
 legstr14 = 'second-order transition point'
 legend(legstr11,legstr12,legstr13,legstr14,'Interpreter','latex', 'Location', 'southeast','FontSize',16)
grid on

subplot(1,2,2)
pltx = linspace(0,1/betajend,100);
%plotdat = inter(pltx);
plotdat = exp(p(2))* pltx.^(p(1));
plot(pltx,plotdat, aaalinear(:,1),aaalinear(:,2),"o","MarkerSize",4)
hold on
plot(aaaremain(:,1),aaaremain(:,2),"-*","MarkerSize",4)
hold on
plot(1/betajend,kappaend,"*","Color","black","MarkerSize",4)
set(gca,'FontSize',15)

secondorderPtstring = append( '(', num2str(1/betajend,'%.3g'), ',', num2str(kappaend,'%.3g'), ')');
text(1/betajend-0.005,kappaend-0.025,secondorderPtstring,'FontSize', 12)
% text(0.5/betajend -0.03,0.5*kappaend+0.25,'p-spin phase','FontSize', 16)
% text(0.5/betajend,0.5*kappaend,'SYK phase','FontSize', 16)
text(0.5/betajend -0.1,0.5*kappaend+0.5,'quasi-integrable phase','FontSize', 18)
text(0.5/betajend,0.5*kappaend+0.1,'chaotic phase','FontSize', 18)

xlim([0 1.2/betajend])
ylim([0 kappaend*1.5])
%titlestr2 = 'Phase diagram, linear scale';
titlestr2 = append('Phase diagram ($l_{nz}$ =', num2str(lnz,'%.4g'), ', $l_{zz}$ = ', num2str(lzz,'%.4g'),'), linear scale');
title(titlestr2,'Interpreter','latex','FontSize',20);
xlabel('$1/(\beta \mathbf{J})$','Interpreter','latex','FontSize',20);
 ylabel('$\kappa$','Interpreter','latex','FontSize',20);
  legstr21 = append('$\kappa =', num2str( exp(p(2)),'%.4f'),'(1/\beta \mathbf{J})^{', num2str(p(1),'%.4f'),'}$');
 legstr22 = 'first-order transition points at low T';
 legstr23 = 'other first-order transition points';
 legstr24 = 'second-order transition point'
  legend(legstr21,legstr22,legstr23,legstr24,'Interpreter','latex','Location', 'southeast','FontSize',16)

  
grid on
%%
% 
% figure;
% pltx = linspace(0,1/betajend,100);
% %plotdat = inter(pltx);
% plotdat = exp(p(2))* pltx.^(p(1));
% plot(pltx,plotdat, aaalinear(:,1),aaalinear(:,2),"o","MarkerSize",4)
% hold on
% plot(aaaremain(:,1),aaaremain(:,2),"-*","MarkerSize",4)
% hold on
% plot(1/betajend,kappaend,"*","Color","black","MarkerSize",4)
% secondorderPtstring = append( '(', num2str(1/betajend,'%.3g'), ',', num2str(kappaend,'%.3g'), ')');
% text(1/betajend-0.0,kappaend-0.02,secondorderPtstring)
% text(0.5/betajend -0.03,0.5*kappaend+0.25,'p-spin phase','FontSize', 16)
% text(0.5/betajend,0.5*kappaend,'SYK phase','FontSize', 16)
% 
% xlim([0 0.01])
% ylim([0 0.2])
% %titlestr2 = 'Phase diagram, linear scale';
% titlestr2 = append('Phase diagram ($l_{nz}$ =', num2str(lnz,'%.4g'), ', $l_{zz}$ = ', num2str(lzz,'%.4g'),'), linear scale');
% title(titlestr2,'Interpreter','latex')
% xlabel('$1/(\beta \mathbf{J})$','Interpreter','latex')
%  ylabel('$\kappa$','Interpreter','latex')
%   legstr21 = append('$\kappa =', num2str( exp(p(2)),'%.4f'),'(1/\beta \mathbf{J})^{', num2str(p(1),'%.4f'),'}$');
%  legstr22 = 'first-order transition points at low T';
%  legstr23 = 'other first-order transition points';
%  legstr24 = 'second-order transition point'
%   legend(legstr21,legstr22,legstr23,legstr24,'Interpreter','latex','Location', 'southeast')
% 
%   set(gca,'FontSize',14)
%% just the linear scale plot

figure;
pltx = linspace(0,1/betajend,100);
%plotdat = inter(pltx);
plotdat = exp(p(2))* pltx.^(p(1));
plot(pltx,plotdat, aaalinear(:,1),aaalinear(:,2),"o","MarkerSize",4)
hold on
plot(aaaremain(:,1),aaaremain(:,2),"-*","MarkerSize",4)
hold on
plot(1/betajend,kappaend,"*","Color","black","MarkerSize",4)
 set(gca,'FontSize',16);

secondorderPtstring = append( '(', num2str(1/betajend,'%.3g'), ',', num2str(kappaend,'%.3g'), ')');
text(1/betajend-0.0,kappaend-0.02,secondorderPtstring,'FontSize', 16)
text(0.5/betajend -0.03,0.5*kappaend+0.25,'quasi-integrable phase','FontSize', 20)
text(0.5/betajend,0.5*kappaend,'chaotic phase','FontSize', 20)
xlim([0 0.1])
ylim([0 0.6])

%titlestr2 = 'Phase diagram, linear scale';
titlestr2 = append('Phase diagram');
title(titlestr2,'FontSize', 20)
xlabel('$1/(\beta \mathbf{J})$','Interpreter','latex')
 ylabel('$\kappa$','Interpreter','latex','FontSize', 20)
  legstr21 = append('$\kappa =', num2str( exp(p(2)),'%.4f'),'(1/\beta \mathbf{J})^{', num2str(p(1),'%.4f'),'}$');
 legstr22 = 'first-order transition points at low T';
 legstr23 = 'other first-order transition points';
 legstr24 = 'second-order transition point'
  legend(legstr21,legstr22,legstr23,legstr24,'Interpreter','LaTeX','Location', 'southeast','FontSize', 20)
grid on
 % set(gca,'FontSize',14)
%%

figure;
pltx = linspace(0,1/betajend,100);

plot(aaa(:,1),aaa(:,2),"-*","Color","#EDB120","MarkerSize",4)
hold on
plot(1/betajend,kappaend,"*","Color","black","MarkerSize",6)
%secondorderPtstring = append( '(', num2str(1/betajend,'%.3g'), ',', num2str(kappaend,'%.3g'), ')');
%text(1/betajend-0.0,kappaend-0.02,secondorderPtstring,'FontSize', 16)
text(0.5/betajend -0.03,0.5*kappaend+0.25,'quasi-integrable phase','FontSize', 20)
text(0.5/betajend,0.5*kappaend,'chaotic phase','FontSize', 20)

xlim([0 0.1])
ylim([0 0.6])
%titlestr2 = 'Phase diagram, linear scale';
titlestr2 = append('Phase diagram');
title(titlestr2,'FontSize', 20)
xlabel('$1/(\beta \mathbf{J})$','Interpreter','latex','FontSize', 20)
 ylabel('$\kappa$','Interpreter','latex','FontSize', 28)
 legstr22 = 'first-order transition points';
 legstr24 = 'second-order transition point'
  legend(legstr22,legstr24,'Interpreter','LaTeX','Location', 'southeast','FontSize', 20)
grid on
%%
% interinverse= @(x) interp1(1./aaa(:,1),aaa(:,2),x,'spline','extrap')
% 
% pltxinverse = linspace(13.5,300,100);
% plotdatinverse = interinverse(pltxinverse);
% figure
% plot(pltx,plotdat, aaa(:,1),aaa(:,2),"o")
%  xlim([0 1/12])
%  ylim([0 .6])
%  xlabel('1/(\beta J)')
%  ylabel('k')
% plot(1./pltxinverse,plotdatinverse, aaa(:,1),aaa(:,2),"o")

%% 

% interinverse= @(x) interp1(1./aaa(:,1),aaa(:,2),x,'spline','extrap')
% 
% pltxinverse = linspace(13.5,160,100);
% plotdatinverse = interinverse(pltxinverse);
% 
% figure
% plot(pltxinverse,plotdatinverse, 1./aaa(:,1),aaa(:,2),"o")
%  xlim([10, 160])
%  ylim([0 1])
%  xlabel('\beta J')
%  ylabel('k')