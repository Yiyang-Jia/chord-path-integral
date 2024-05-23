
%clear all
global betaj k l lz v u 

l=1.;
lz=0;
betaj=1500;
% note for decreasing k solutions, initial mesh size can be important
trialpoint =0.05; %trial starting point for searching 1st order transition value of k
decreasingkappamesh=linspace(-0.5, 0.5, 91);% note for decreasing k solutions, initial mesh size can be important

k=0.;

syms f(x) %use this to generate initial guess for g_n at k=0.
f(x) = pi*x/cos(pi* x/2) - betaj;
v = eval(vpasolve(f==0, x, [0 1]));



 options = bvpset('RelTol',1e-4, 'Nmax', 5000);
xmesh=linspace(-0.5, 0.5, 103);
%solinitk0 = bvpinit(xmesh, [-30; 0;-4;0]);  %this is a better guess for lz ~= 0 and lz<<1
solinitk0 = bvpinit(xmesh, @exactk0);
solk0= bvp4c(@bvpfun, @bcfun, solinitk0,options);
%% 

% plot(solk0.x,solk0.y(1,:),solk0.x,solk0.y(3,:))
%  action(solk0)/sykaction(betaj)  %->ratio is always 2 for different
%  temperatures, seems to be an overall factor missing

actionlistincreasek = [k, action(solk0)];  %values of action, with k increasing from 0 
soliterateInc = solk0;


npoints = 600;

for n = 1: npoints
         lastwarn('') % Clear last warning message
    k =  k+ 1/npoints; 
    % if k >0.55
    %      options = bvpset('RelTol',1e-5, 'Nmax', 20000);
    % end
  try
    solintk0 =bvpinit(soliterateInc.x, @(u) interp1(soliterateInc.x, transpose(soliterateInc.y), u));%https://www.mathworks.com/matlabcentral/answers/142088-initial-guess-for-bvp4c
    soliterateInc= bvp4c(@bvpfun, @bcfun, solinitk0,options);
    actionlistincreasek = [actionlistincreasek; k action(soliterateInc)]
    [warnMsg, warnId] = lastwarn;
    warnIncreasek="";
    errIncreasek = "";
     if ~isempty(warnMsg)%break the loop is there's a warning message (poor convergence when solving ode)
     warnIncreasek = "Warning: increasing k sequence, k= "+ k + ". " + warnMsg;
       break 
     end
  catch err
       errIncreasek = "Error: increasing k sequence, k= "+ k + ". " + err.message;
      break
  end 
end
actionlistincreasek(end,:) = []; %remove the last result when the loop breaks, this is a poorly converged result

%% 


k=1;
%decreasingkappamesh=linspace(-0.5, 0.5, 121);% note for decreasing k solutions, initial mesh size can be important
xmesh=decreasingkappamesh;
%npoints = 100;

%options = bvpset('RelTol',1e-4, 'Nmax', 10000);

if lz ~= 0%use this to generate initial guess for g_z at k=1, if lz is not zero
  syms g(x)
  g(x) = pi*x/cos(pi* x/2) - sqrt(lz)*betaj;
  u = eval(vpasolve(g==0, x, [0 1]));  
  %u = 0.99999999; %maybe use this if betaj is extremely large
 solinitk1 = bvpinit(xmesh, @exactk1lz);
else
 solinitk1 = bvpinit(xmesh, @exactk1);
end

%  solinitk1 = bvpinit(xmesh, @exactk1);
%solinitk1 = bvpinit(xmesh, [-0.01;0;-35;0]); %Actually this is better guess for lz ~=0 and <<1 
solk1 = bvp4c(@bvpfun, @bcfun, solinitk1, options);
actionlistdecreasek = [k, action(solk1)];  %actions, with k decreasing from 1

soliterateDec =solk1;
%%
for n = 1: npoints
       lastwarn('') % Clear last warning message
      k =  k - 1/npoints; 
      if k<0
          break
      end
        
    %   if k <0.56
    %      options = bvpset('RelTol',1e-5, 'Nmax', 20000);
    % end
      try
      solintk1 =bvpinit(soliterateDec.x, @(u) interp1(soliterateDec.x, transpose(soliterateDec.y), u));%https://www.mathworks.com/matlabcentral/answers/142088-initial-guess-for-bvp4c
      soliterateDec= bvp4c(@bvpfun, @bcfun, solinitk1,options);
      actionlistdecreasek = [actionlistdecreasek; k action(soliterateDec)]
    
       [warnMsg, warnId] = lastwarn;
       warnDecreasek = "";
       errDecreasek = "" ;
      if ~isempty(warnMsg)%break the loop is there's a warning message (poor convergence when solving ode)
     warnDecreasek = "Warning: decreasing k sequence, k= "+ k + ". " + warnMsg;
       break 
      end
  catch err
      errDecreasek = "Error: decreasing k sequence, k= "+ k + ". " + err.message;
      break
   end

end
actionlistdecreasek(end,:) = [];
 
 fprintf( '%s\n',warnIncreasek)
   fprintf('%s\n', errIncreasek)

 fprintf('%s\n',warnDecreasek)
 fprintf('%s\n', errDecreasek)
% 
%% Plot the last solutions before numerics break down
% 
% figure;
%  subplot(2,1,1)
%    plot(soliterateDec.x, soliterateDec.y(1,:),soliterateDec.x, soliterateDec.y (3,:))
%    title('last point of decreasing kappa')
%    legend('gn','gz')
%     subplot(2,1,2)
%   plot(soliterateInc.x, soliterateInc.y(1,:),soliterateInc.x, soliterateInc.y (3,:))
%      title('last point of increasing kappa')
%         legend('gn','gz')

%% find subdominant actions  
%    %we roughly know the subdominant solution at  betaj=65 and kappa=0.24
% k=0.24;
% 
%  options = bvpset('RelTol',1e-5, 'Nmax', 5000);
% xmesh=linspace(-0.5, 0.5, 103);
% solinitk0 = bvpinit(xmesh, [-1; 0;-6;0]);    %we roughly know the subdominant solution at  betaj=65 and kappa=0.24
% solk0= bvp4c(@bvpfun, @bcfun, solinitk0,options);
% actionlistincreaseksubdom = [k, action(solk0)];  %values of action, with k increasing from 0 
% soliterateInc = solk0;
% 
% 
% npoints = 200;
% 
% for n = 1: npoints
%          lastwarn('') % Clear last warning message
%     k =  k+ 1/npoints; 
%     % if k >0.55
%     %      options = bvpset('RelTol',1e-5, 'Nmax', 20000);
%     % end
%   try
%     solintk0 =bvpinit(soliterateInc.x, @(u) interp1(soliterateInc.x, transpose(soliterateInc.y), u));%https://www.mathworks.com/matlabcentral/answers/142088-initial-guess-for-bvp4c
%     soliterateInc= bvp4c(@bvpfun, @bcfun, solinitk0,options);
%     actionlistincreaseksubdom = [actionlistincreaseksubdom; k action(soliterateInc)]
%     [warnMsg, warnId] = lastwarn;
%     warnIncreasek="";
%     errIncreasek = "";
%      if ~isempty(warnMsg)%break the loop is there's a warning message (poor convergence when solving ode)
%      warnIncreasek = "Warning: increasing k sequence, k= "+ k + ". " + warnMsg;
%        break 
%      end
%   catch err
%        errIncreasek = "Error: increasing k sequence, k= "+ k + ". " + err.message;
%       break
%   end 
% end
% actionlistincreaseksubdom(end,:) = []; %remove the last result when the loop breaks, this is a poorly converged result
% 
% 
% k=0.24-1/npoints;
% 
% solinitk1 = bvpinit(xmesh, [-1; 0;-6;0]);    %we roughly know the subdominant solution at  betaj=65 and kappa=0.24
% solk1 = bvp4c(@bvpfun, @bcfun, solinitk1, options);
% actionlistdecreaseksubdom = [k, action(solk1)];  %actions, with k decreasing from 1
% 
% soliterateDec =solk1;
% for n = 1: npoints
%        lastwarn('') % Clear last warning message
%       k =  k - 1/npoints; 
%       if k<0
%           break
%       end
% 
%     %   if k <0.56
%     %      options = bvpset('RelTol',1e-5, 'Nmax', 20000);
%     % end
%       try
%       solintk1 =bvpinit(soliterateDec.x, @(u) interp1(soliterateDec.x, transpose(soliterateDec.y), u));%https://www.mathworks.com/matlabcentral/answers/142088-initial-guess-for-bvp4c
%       soliterateDec= bvp4c(@bvpfun, @bcfun, solinitk1,options);
%       actionlistdecreaseksubdom = [actionlistdecreaseksubdom; k action(soliterateDec)]
% 
%        [warnMsg, warnId] = lastwarn;
%        warnDecreasek = "";
%        errDecreasek = "" ;
%       if ~isempty(warnMsg)%break the loop is there's a warning message (poor convergence when solving ode)
%      warnDecreasek = "Warning: decreasing k sequence, k= "+ k + ". " + warnMsg;
%        break 
%       end
%   catch err
%       errDecreasek = "Error: decreasing k sequence, k= "+ k + ". " + err.message;
%       break
%    end
% 
% end
% actionlistdecreaseksubdom(end,:) = [];


%%
   % plot subdominant actions betaj=65
% actionlistsubdom = sortrows([actionlistdecreaseksubdom;actionlistincreaseksubdom])
% figure;
% subplot(1,2,1)
%  plot(actionlistincreasek(:,1),actionlistincreasek(:,2),actionlistdecreasek(:,1),actionlistdecreasek(:,2),'LineWidth',1)
% hold on
%  plot(actionlistsubdom(:,1),actionlistsubdom(:,2),'Color',"#7E2F8E",'LineWidth',1)
%   %titlestr= append("Actions for varying $\kappa$ ($\beta \mathbf{J}$ = " , num2str(betaj, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')+")") 
%    titlestr= append("Actions for varying $\kappa$ ($\beta \mathbf{J}$ = " , num2str(betaj, '%.4g'), ")") 
% 
%   title(titlestr,'Interpreter','latex',FontSize=18)
%   xlabel('$\kappa$','Interpreter','latex',FontSize=20)
%   ylabel('Action', 'Interpreter','latex', FontSize=16)
%  lgstr1 ='chaotic phase';
%   lgstr2 = 'integrable phase';
%    lgstr3='subdominant saddle';
%     legend(lgstr1,lgstr2,lgstr3,'Interpreter','latex','Location', 'southwest',FontSize=18)
% 
% subplot(1,2,2)
%  plot(actionlistincreasek(:,1),actionlistincreasek(:,2),actionlistdecreasek(:,1),actionlistdecreasek(:,2),'LineWidth',1.5)
% hold on
%  plot(actionlistsubdom(:,1),actionlistsubdom(:,2),'Color',"#7E2F8E",'LineWidth',1.5)
% 
%  xlim([0.1,0.4]);
% 
%  %titlestr= append("Actions for varying $\kappa$ ($\beta \mathbf{J}$ = " , num2str(betaj, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')+")") 
%    titlestr= append("Actions for varying $\kappa$ ($\beta \mathbf{J}$ = " , num2str(betaj, '%.4g'), "), zoomed-in") 
% 
%   title(titlestr,'Interpreter','latex',FontSize=18)
%   xlabel('$\kappa$','Interpreter','latex',FontSize=20)
%   ylabel('Action', 'Interpreter','latex', FontSize=16)
%  lgstr1 ='chaotic phase';
%   lgstr2 = 'integrable phase';
%    lgstr3='subdominant saddle';
%     legend(lgstr1,lgstr2,lgstr3,'Interpreter','latex','Location', 'southwest',FontSize=18)

%%
 %   %write to file for betaj=65
 %      filename1 = append('actionsbetaj=65chaotic.txt');
 % dlmwrite(filename1 ,actionlistincreasek,'-append'); 
 %   filename2 = append('actionsbetaj=65integrable.txt');
 % dlmwrite(filename2 ,actionlistdecreasek,'-append'); 
 %   filename3 = append('actionsbetaj=65subdom.txt');
 % dlmwrite(filename3 ,actionlistsubdom,'-append'); 
  %%
% figure;
%  plot(actionlistincreasek(:,1),actionlistincreasek(:,2),actionlistdecreasek(:,1),actionlistdecreasek(:,2),'LineWidth',1)
%  titlestr= append("Action as a function of $\kappa$ ($\beta \mathbf{J}$ = " , num2str(betaj, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')+")") 
%  title(titlestr,'Interpreter','latex',FontSize=14)
%   xlabel('$\kappa$','Interpreter','latex',FontSize=20)
%   ylabel('Action', 'Interpreter','latex', FontSize=16)

% truncate curves to [0,0.5] and [0.5,1]
% figure;
%  plot(actionlistincreasek(1:npoints/2+4,1),actionlistincreasek(1:npoints/2+4,2),actionlistdecreasek(1:npoints/2,1),actionlistdecreasek(1:1:npoints/2,2),'LineWidth',1)
%  titlestr= append("Action as a function of $\kappa$ ($\beta \mathbf{J}$ = " , num2str(betaj, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')+")") 
%  title(titlestr,'Interpreter','latex',FontSize=14)
%   xlabel('$\kappa$','Interpreter','latex',FontSize=20)
%   ylabel('Action', 'Interpreter','latex', FontSize=16)
% 
%%
% 
% kappamesh =linspace(0,1,300);
% betaj100=100;
%  figure;
% subplot(1,3,1)
% plot(actionlistincreasek100(:,1),actionlistincreasek100(:,2),'.','LineWidth',0.5,'Color',"#0072BD")
% hold on
% plot(kappamesh, -2*betaj100* sqrt(1-kappamesh.^2),'LineWidth',0.5,'Color',"#0072BD")
% hold on
% plot(actionlistdecreasek100(:,1),actionlistdecreasek100(:,2),'.','LineWidth',0.5,'Color',"#D95319")
% hold on
% plot(kappamesh,-1/2*(betaj100*kappamesh).^2,'LineWidth',0.5,'Color',"#D95319")
%  %titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')) 
%   titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g')) 
% 
%  title(titlestr,'Interpreter','latex',FontSize=20)
%   xlabel('$\kappa$','Interpreter','latex',FontSize=20)
%   ylabel('Action', 'Interpreter','latex', FontSize=20)
% lgstr1 ='chaotic phase';
%   lgstr3='integrable phase';
%   lgstr2 = 'chaotic phase approx.';
%   lgstr4= 'integrable phase approx.';
%   legend(lgstr1,lgstr2,lgstr3,lgstr4,'Interpreter','latex','Location', 'southwest',FontSize=18)
%     % set(gca,'FontSize',18)
% 
% subplot(1,3,2)
% plot(actionlistincreasek100(:,1),actionlistincreasek100(:,2),'.','LineWidth',0.5,'Color',"#0072BD")
% hold on
% plot(kappamesh, -2*betaj100* sqrt(1-kappamesh.^2),'LineWidth',1,'Color',"#0072BD")
% hold on
% plot(actionlistdecreasek100(:,1),actionlistdecreasek100(:,2),'.','LineWidth',0.5,'Color',"#D95319")
% hold on
% plot(kappamesh,-1/2*(betaj100*kappamesh).^2,'LineWidth',1,'Color',"#D95319")
%  %titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')) 
%   titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g'),', zoomed-in') 
% 
%  title(titlestr,'Interpreter','latex',FontSize=20)
%   xlabel('$\kappa$','Interpreter','latex',FontSize=20)
%   ylabel('Action', 'Interpreter','latex', FontSize=20)
% lgstr1 ='chaotic phase';
%   lgstr3='integrable phase';
%   lgstr2 = 'chaotic phase approx.';
%   lgstr4= 'integrable phase approx.';
%   legend(lgstr1,lgstr2,lgstr3,lgstr4,'Interpreter','latex','Location', 'southwest',FontSize=18)
%   xlim([0.15,0.25])
% 
%     % set(gca,'FontSize',18)
% 
% subplot(1,3,3)
% plot(actionlistincreasek(:,1),actionlistincreasek(:,2),'.','LineWidth',0.5,'Color',"#0072BD")
% hold on
% plot(kappamesh, -2*betaj* sqrt(1-kappamesh.^2),'LineWidth',1,'Color',"#0072BD")
% hold on
% plot(actionlistdecreasek(:,1),actionlistdecreasek(:,2),'.','LineWidth',0.5,'Color',"#D95319")
% hold on
% plot(kappamesh,-1/2*(betaj*kappamesh).^2,'LineWidth',1,'Color',"#D95319")
%  %titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')) 
%   titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj, '%.4g'),', zoomed-in') 
% 
%  title(titlestr,'Interpreter','latex',FontSize=20)
%   xlabel('$\kappa$','Interpreter','latex',FontSize=20)
%   ylabel('Action', 'Interpreter','latex', FontSize=20)
% lgstr1 ='chaotic phase';
%   lgstr3='integrable phase';
%   lgstr2 = 'chaotic phase approx.';
%   lgstr4= 'integrable phase approx.';
%   legend(lgstr1,lgstr2,lgstr3,lgstr4,'Interpreter','latex','Location', 'southwest',FontSize=18)
%   xlim([0,0.1])
    % set(gca,'FontSize',18)
 % 
 % actionlistdecreasek1000 =  actionlistdecreasek;
 % actionlistincreasek1000=actionlistincreasek;

 %% action curve with Joseif's correction

kappamesh =linspace(0,1,300);
betaj100=100;
 figure;
subplot(1,3,1)
plot(actionlistincreasek100(:,1),actionlistincreasek100(:,2),'.','LineWidth',0.5,'Color',"#0072BD")
hold on
plot(kappamesh,0*kappamesh-2*betaj100,'LineWidth',0.5,'Color',"#0072BD")
hold on
plot(actionlistdecreasek100(:,1),actionlistdecreasek100(:,2),'.','LineWidth',0.5,'Color',"#D95319")
hold on
plot(kappamesh,-1/2*(betaj100*kappamesh).^2,'LineWidth',0.5,'Color',"#D95319")
 %titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')) 
 
 set(gca,'FontSize',16)

 titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g')) 

 title(titlestr,'Interpreter','latex',FontSize=20)
  xlabel('$\kappa$','Interpreter','latex',FontSize=20)
  ylabel('Action', 'Interpreter','latex', FontSize=20)
lgstr1 ='chaotic phase';
  lgstr3='integrable phase';
  lgstr2 = 'chaotic phase approx.';
  lgstr4= 'integrable phase approx.';
  legend(lgstr1,lgstr2,lgstr3,lgstr4,'Interpreter','latex','Location', 'southwest',FontSize=18)
    % set(gca,'FontSize',18)
    grid on

subplot(1,3,2)
plot(actionlistincreasek100(:,1),actionlistincreasek100(:,2),'.','LineWidth',0.5,'Color',"#0072BD")
hold on
plot(kappamesh, 0*kappamesh-2*betaj100,'LineWidth',1,'Color',"#0072BD")
hold on
plot(actionlistdecreasek100(:,1),actionlistdecreasek100(:,2),'.','LineWidth',0.5,'Color',"#D95319")
hold on
plot(kappamesh,-1/2*(betaj100*kappamesh).^2,'LineWidth',1,'Color',"#D95319")
 %titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')) 
  
 set(gca,'FontSize',14)

 titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g'),', zoomed-in') 

 title(titlestr,'Interpreter','latex',FontSize=20)
  xlabel('$\kappa$','Interpreter','latex',FontSize=20)
  ylabel('Action', 'Interpreter','latex', FontSize=20)
lgstr1 ='chaotic phase';
  lgstr3='integrable phase';
  lgstr2 = 'chaotic phase approx.';
  lgstr4= 'integrable phase approx.';
  legend(lgstr1,lgstr2,lgstr3,lgstr4,'Interpreter','latex','Location', 'southwest',FontSize=18)
  xlim([0.15,0.25])
  ylim([-400,-150])
    grid on


subplot(1,3,3)
plot(actionlistincreasek(:,1),actionlistincreasek(:,2),'.','LineWidth',0.5,'Color',"#0072BD")
hold on
plot(kappamesh, 0*kappamesh-2*betaj,'LineWidth',1,'Color',"#0072BD")
hold on
plot(actionlistdecreasek(:,1),actionlistdecreasek(:,2),'.','LineWidth',0.5,'Color',"#D95319")
hold on
plot(kappamesh,-1/2*(betaj*kappamesh).^2,'LineWidth',1,'Color',"#D95319")

  
set(gca,'FontSize',14)
%titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')) 
  titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj, '%.4g'),', zoomed-in') 
 title(titlestr,'Interpreter','latex',FontSize=20)
  xlabel('$\kappa$','Interpreter','latex',FontSize=20)
  ylabel('Action', 'Interpreter','latex', FontSize=20)
lgstr1 ='chaotic phase';
  lgstr3='integrable phase';
  lgstr2 = 'chaotic phase approx.';
  lgstr4= 'integrable phase approx.';
  legend(lgstr1,lgstr2,lgstr3,lgstr4,'Interpreter','latex','Location', 'southwest',FontSize=18)
  xlim([0,0.1])
    grid on
 %% plot just the betaj=100 case

kappamesh =linspace(0,1,300);
betaj100=100;
 figure;
plot(actionlistincreasek100sparse(:,1),actionlistincreasek100sparse(:,2),'o','MarkerSize',5,'Color',"#0072BD")
hold on
plot(kappamesh, 0*kappamesh-2*betaj100,'-','LineWidth',1,'Color',"#0072BD")
hold on
plot(actionlistdecreasek100sparse(:,1),actionlistdecreasek100sparse(:,2),'o','MarkerSize',5,'Color',"#D95319")
hold on
plot(kappamesh,-1/2*(betaj100*kappamesh).^2,'-','LineWidth',0.5,'Color',"#D95319")

    set(gca,'FontSize',14)

%titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')) 
  titlestr= append("Action for varying $\kappa$, $\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g')) 

 title(titlestr,'Interpreter','latex',FontSize=20)
  xlabel('$\kappa$','Interpreter','latex',FontSize=20)
  ylabel('Action', 'Interpreter','latex', FontSize=20)
lgstr1 ='chaotic phase numerical';
  lgstr3='quasi-integrable phase numerical';
  lgstr2 = 'chaotic phase approx.';
  lgstr4= 'quasi-integrable phase approx.';
  legend(lgstr1,lgstr2,lgstr3,lgstr4,'Interpreter','latex','Location', 'southwest',FontSize=20)
  grid on
%%
% kappamesh =linspace(0,1,300);
%  figure;
% plot(actionlistincreasek(:,1),actionlistincreasek(:,2),'o','MarkerSize',5,'Color',"#0072BD")
% hold on
% plot(kappamesh, 0*kappamesh-2*betaj,'-','LineWidth',1,'Color',"#0072BD")
% hold on
% plot(actionlistdecreasek(:,1),actionlistdecreasek(:,2),'o','MarkerSize',5,'Color',"#D95319")
% hold on
% plot(kappamesh,-1/2*(betaj*kappamesh).^2,'-','LineWidth',0.5,'Color',"#D95319")
%  %titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj100, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')) 
%   titlestr= append("Action for varying $\kappa$, $\beta \mathbf{J}$ = " , num2str(betaj, '%.4g')) 
% 
%  title(titlestr,'Interpreter','latex',FontSize=20)
%   xlabel('$\kappa$','Interpreter','latex',FontSize=20)
%   ylabel('Action', 'Interpreter','latex', FontSize=20)
% lgstr1 ='chaotic phase numerical';
%   lgstr3='integrable phase numerical';
%   lgstr2 = 'chaotic phase approx.';
%   lgstr4= 'integrable phase approx.';
%   legend(lgstr1,lgstr2,lgstr3,lgstr4,'Interpreter','latex','Location', 'southwest',FontSize=20)
%       grid on

    %% 
% 
% 
% %trialpoint =0.43; %find this trial point by eye
%  transitionPts = findcritk(actionlistincreasek,actionlistdecreasek ,trialpoint)
% 
% 
%  lstr =strsplit(num2str(l),'.');
%  lzstr = strsplit(num2str(lz),'.');
%  if length(lstr) > 1
%      lstrcom = lstr{1} + "p" + lstr{2};
%  else
%      lstrcom = lstr{1};
%  end
% 
%  if length(lzstr) > 1
%      lzstrcom = lzstr{1} + "p" + lzstr{2};
%  else
%      lzstrcom = lzstr{1};
%  end
%  filename = append('tranptsHigherAccuracyl=', lstrcom,"lz=", lzstrcom,".txt");
%  dlmwrite( filename ,transitionPts,'-append');  %in the format of [1/betaj  kc]


%% 

  [p,S]=polyfit(actionlistdecreasek(:,1),actionlistdecreasek(:,2),1);

%%




function dydx = bvpfun(x,y) % equation being solved
global betaj k l lz
dydx=zeros(4,1);
dydx = [y(2)
        2*betaj^2*(1-k^2)*exp(y(1)+l*y(3))
        y(4)
         2*betaj^2* k^2*exp(l*y(1) + lz* y(3))];
end

function res = bcfun(ya,yb) % boundary conditions 
res = [ya(1)
       yb(1)
       ya(3)
       yb(3)];
end

function ggg = exactk0(z) % exact solution for gn when k=0
global betaj v
ggg = [2 * log(cos(pi*v/2)/cos(pi*v*z))
       2*pi*v * tan(pi*v*z)
       0
       0];
end

function ggg = exactk1lz(z) % exact solution for gz when k=1 and lz ~= 0
global betaj u lz
ggg = [2/lz * log(cos(pi*u/2)/cos(pi*u*z))
       2*pi*u/lz * tan(pi*u*z)
       0
       0];
end

function g = exactk1(x) % exact solution when k=1
global betaj 
g = [0
    0
    betaj^2 * x^2 - betaj^2/4
    2 *betaj^2 *x];
end



function s = action(sol)  %evaluate action
global betaj k l lz 
x = sol.x;
x = x(:);
gn = sol.y(1,:);
gn = gn(:);
gz = sol.y(3,:);
gz = gz(:);
integrand =0.5.* betaj^2.*(0.5 - x).* ((1-k^2).*exp(gn+l*gz).* (gn-2) +l*(1-k^2).*exp(gn+l*gz).* gz + k^2.* exp(l*gn+lz*gz) .* (l*gn-2) + lz *k^2.* exp(l*gn+lz*gz) .* gz);
   integrandinterp =  @(u) interp1(x, integrand, u);
    s = integral(integrandinterp, -0.5, 0.5);
%     fplot(integrandinterp,[-0.5,0.5])
end

function ssyk = sykaction(betaj) %action of pure syk model
global v
% eval(v - pi)
ssyk = eval( pi^2*v^2/4 - pi * v* tan(pi * v /2));
end

function crit = findcritk(kaction1, kaction2,trialpoint) %action of pure syk model
global betaj 
curve1 = @(u) interp1(kaction1(:,1), kaction1(:,2), u,'spline','extrap');
curve2 = @(u) interp1(kaction2(:,1), kaction2(:,2), u,'spline','extrap');
diff = @(x) curve1(x) - curve2(x);
kc = fsolve(diff, trialpoint);
crit = [1/betaj kc];
end

