
aaa = readmatrix('tranptsHigherAccuracyl=1lz=0p005.txt');

lnz=1;
lzz=0.005;
kappazero = sqrt(lzz/(1+lzz));
slopezero =pi^2/4 *(1/lzz-1)/(1/sqrt(lzz)+sqrt(lzz));
betajbegin=1/aaa(1,1);
betajend = 1/aaa(end,1);
kappaend = aaa(end, 2);
datasplitpoint =10; %only low temperature point obey linear relation
aaalinear =aaa(1:datasplitpoint,:);
aaaremain = aaa(datasplitpoint+1:end,:);


%inter= @(x) interp1(aaa(:,1),aaa(:,2),x,'spline','extrap');

% template = hgload('phasediagramloglog.fig');
% ax=gca;
% ha1 = findobj(template,'Type','Legend');
% 
% 
   % [p,S, mu]=polyfit(aaalinear(:,1),aaalinear(:,2),1);
   [p,S]=polyfit(aaalinear(:,1),aaalinear(:,2),1);
  pltx2 = linspace(0,1.5*aaalinear(end,1),100);
%     [y,delta] = polyval(p,pltx2,S,mu);
% plot(pltx2,y+2*delta,'m--',pltx2,y-2*delta,'m--')
% hold on
% plot(aaalinear(:,1),aaalinear(:,2),"o","MarkerSize",4)
 %%
 ff = figure;
  
subplot(1,2,1)
pltx = linspace(0,1/betajend,100);
%plotdat = inter(pltx);
%plotdat = exp(p(2))* pltx.^(p(1));
%plot(pltx,plotdat, aaalinear(:,1),aaalinear(:,2),"o","MarkerSize",4)
plot(aaa(1:end-1,1),aaa(1:end-1,2),"-o",'Color','Red',"MarkerSize",4)
% hold on
% plot(aaaremain(:,1),aaaremain(:,2),"-*","MarkerSize",4)
hold on
plot(1/betajend,kappaend,"*","Color","black","MarkerSize",4)
hold on
plot(aaa(end-1:end,1),aaa(end-1:end,2),"-",'Color','Red',"MarkerSize",4)
  set(gca,'FontSize',16)

secondorderPtstring = append( '(', num2str(1/betajend,'%.3g'), ',', num2str(kappaend,'%.3g'), ')');
text(1/betajend-0.004,kappaend-0.03,secondorderPtstring,'FontSize', 11)
text(0.5/betajend -0.03,0.5*kappaend+0.25,'polarized phase','FontSize', 16)
text(0.5/betajend,0.5*kappaend,'chaotic phase','FontSize', 16)

xlim([0 1.2/betajend])
ylim([0 kappaend*1.5])
%titlestr2 = 'Phase diagram, linear scale';
titlestr2 = append('Phase diagram ($l_{nz}$ =', num2str(lnz,'%.4g'), ', $l_{zz}$ = ', num2str(lzz,'%.4g'),')');
title(titlestr2,'Interpreter','latex')
xlabel('$1/(\beta \mathbf{J})$','Interpreter','latex')
 ylabel('$\kappa$','Interpreter','latex')
%  legstr21 = append('$\kappa =', num2str( exp(p(2)),'%.4f'),'(1/\beta \mathbf{J})^{', num2str(p(1),'%.4f'),'}$');
 legstr22 = 'first-order transition points';
% legstr23 = 'other first-order transition points';
 legstr24 = 'second-order transition point'
  legend(legstr22,legstr24,'Interpreter','latex','Location', 'southeast','FontSize', 16)
grid on

subplot(1,2,2)
plot(pltx2, p(1)*pltx2+ p(2),'-','Color','Blue');
hold on
plot(aaalinear(:,1),aaalinear(:,2),"o",'Color','Red',"MarkerSize",4)
 hold on
 plot(aaaremain(:,1),aaaremain(:,2),"o","MarkerSize",4)
 set(gca,'FontSize',15)
 % hold on
% plot(log(1/betajend),log(kappaend),"*","Color","black","MarkerSize",4)
xlim([0 1.5*aaalinear(end,1)])
ylim([0.03 1.5*aaalinear(end,2)])
%titlestr1 = 'Phase diagram, log-log scale';
titlestr1 = append('Phase diagram ($l_{nz}$ =', num2str(lnz,'%.4g'), ', $l_{zz}$ = ', num2str(lzz,'%.4g'), '), low temperatures');
title(titlestr1,'Interpreter','latex')

text(0.0002,0.08,'polarized phase','FontSize', 16)
text(0.0006,0.06,'chaotic phase','FontSize', 16)

 xlabel('$1/(\beta \mathbf{J} )$','Interpreter','latex')
 ylabel('$\kappa$','Interpreter','latex')
 legstr11 = append('$\kappa=',num2str(p(1),'%.4f'),' (1/\beta \mathbf{J})$ + ', num2str(p(2),'%.4f'));
 legstr12 = 'first-order transition points at low T';
 legstr13 ='other first-order transition points'
 %legstr14 = 'second-order transition point'
 legend(legstr11,legstr12,legstr13,'Interpreter','latex', 'Location', 'southeast','FontSize', 16)
grid on 

%%

figure;
subplot(1,2,1)
plot(log(pltx2), log(p(1)*pltx2+ p(2)),'-','Color','Blue');
hold on
plot(log(aaalinear(:,1)),log(aaalinear(:,2)),"o",'Color','Red',"MarkerSize",4)
 hold on
 plot(log(aaaremain(:,1)),log(aaaremain(:,2)),"o","MarkerSize",4)
% hold on
% plot(log(1/betajend),log(kappaend),"*","Color","black","MarkerSize",4)
xlim([1.5*log(aaalinear(end,1)), -3])
ylim([1.5*log(aaalinear(end,2)),-0])
%titlestr1 = 'Phase diagram, log-log scale';
titlestr1 = append('Phase diagram ($l_{nz}$ =', num2str(lnz,'%.4g'), ', $l_{zz}$ = ', num2str(lzz,'%.4g'), '), low temperatures');
title(titlestr1,'Interpreter','latex')

%text(0.2,0.08,'non-isotropic SYK phase','FontSize', 16)
%text(0.6,0.05,'SYK phase','FontSize', 16)

 xlabel('$\log 1/(\beta \mathbf{J} )$','Interpreter','latex')
 ylabel('$\log \kappa$','Interpreter','latex')
 legstr11 = append('$\log \kappa= \log [',num2str(p(1),'%.4f'),' (1/\beta \mathbf{J})$ + ', num2str(p(2),'%.4f'),']');
 legstr12 = 'first-order transition points at low T';
 legstr13 ='other first-order transition points'
 %legstr14 = 'second-order transition point'
 legend(legstr11,legstr12,legstr13,'Interpreter','latex', 'Location', 'southeast')
 set(gca,'FontSize',14)

 subplot(1,2,2)
plot(log(pltx2), log(p(1)*pltx2),'-','Color','Blue');
hold on
plot(log(aaalinear(:,1)),log(aaalinear(:,2)-p(2)),"o",'Color','Red',"MarkerSize",4)
 hold on
 plot(log(aaaremain(:,1)),log(aaaremain(:,2)-p(2)),"o","MarkerSize",4)
% hold on
% plot(log(1/betajend),log(kappaend),"*","Color","black","MarkerSize",4)
xlim([1.2*log(aaalinear(end,1)), -4])
ylim([2*log(aaalinear(end,2)),0])
%titlestr1 = 'Phase diagram, log-log scale';
titlestr1 = append('Phase diagram ($l_{nz}$ =', num2str(lnz,'%.4g'), ', $l_{zz}$ = ', num2str(lzz,'%.4g'), '), low temperatures');
title(titlestr1,'Interpreter','latex')

text(0.2,0.08,'non-isotropic SYK phase','FontSize', 16)
text(0.6,0.05,'SYK phase','FontSize', 16)

 xlabel('$\log 1/(\beta \mathbf{J} )$','Interpreter','latex')
 ylabel('$\log (\kappa - \kappa_{0,fit})$','Interpreter','latex')
 legstr11 = append('$\log (\kappa - ', num2str(p(2),'%.4f'),' ) =\log [',num2str(p(1),'%.4f'),' (1/\beta \mathbf{J})$ ]');
 legstr12 = 'first-order transition points at low T';
 legstr13 ='other first-order transition points'
 %legstr14 = 'second-order transition point'
 legend(legstr11,legstr12,legstr13,'Interpreter','latex', 'Location', 'southeast')
 set(gca,'FontSize',14)
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