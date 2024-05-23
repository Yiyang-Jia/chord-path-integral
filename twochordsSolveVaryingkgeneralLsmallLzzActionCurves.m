
%clear all
global betaj k l lz v u 

l=1.; %next try l=2, lz=4, this is the p-2p model
lz=0.005;
betaj=500;
% note for decreasing k solutions, initial mesh size can be important
trialpoint =0.05; %trial starting point for searching 1st order transition value of k
decreasingkappamesh=linspace(-0.5, 0.5, 111);% note for decreasing k solutions, initial mesh size can be important

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


npoints = 100;

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
npoints = 100;

%options = bvpset('RelTol',1e-4, 'Nmax', 10000);

% if lz ~= 0%use this to generate initial guess for g_z at k=1, if lz is not zero
%   syms g(x)
%   g(x) = pi*x/cos(pi* x/2) - sqrt(lz)*betaj;
%   u = eval(vpasolve(g==0, x, [0 1]));  
%   %u = 0.99999999; %maybe use this if betaj is extremely large
%  solinitk1 = bvpinit(xmesh, @exactk1lz);
% else
%  solinitk1 = bvpinit(xmesh, @exactk1);
% end

%  solinitk1 = bvpinit(xmesh, @exactk1);
solinitk1 = bvpinit(xmesh, [-0.01;0;-35;0]); %Actually this is better guess for lz ~=0 and <<1 
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
% % 
% figure;
%  subplot(2,1,1)
%    plot(soliterateDec.x, soliterateDec.y(1,:),soliterateDec.x, soliterateDec.y (3,:))
%    title('last point of decreasing kappa')
%    legend('gn','gz')
%     subplot(2,1,2)
%   plot(soliterateInc.x, soliterateInc.y(1,:),soliterateInc.x, soliterateInc.y (3,:))
%      title('last point of increasing kappa')
%         legend('gn','gz')
% 

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
kappamesh =linspace(0,1,100);
betaj1000=1000;
figure;
subplot(1,2,2)
plot(actionlistincreasek1000(:,1),actionlistincreasek1000(:,2),'o',actionlistdecreasek1000(:,1),actionlistdecreasek1000(:,2),'o','MarkerSize',4)
hold on
plot(kappamesh, 0*kappamesh-2*betaj1000+pi^2/2, kappamesh,-2*betaj1000/sqrt(lz)*kappamesh+pi^2/2/lz,'LineWidth',1)
    set(gca,'FontSize',16)

titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj1000, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')) 
 title(titlestr,'Interpreter','latex',FontSize=20)
  xlabel('$\kappa$','Interpreter','latex',FontSize=20)
  ylabel('Action', 'Interpreter','latex', FontSize=20)
lgstr1 ='increasing $\kappa$ numerical';
  lgstr2='decreasing $\kappa$ numerical';
  lgstr3 = 'increasing $\kappa$ analytic approx.';
  lgstr4= 'decreasing $\kappa$ analytic approx.';
  legend(lgstr1,lgstr2,lgstr3,lgstr4,'Interpreter','latex','Location', 'southwest',FontSize=18)
grid on

subplot(1,2,1)
 plot(actionlistincreasek(:,1),actionlistincreasek(:,2),'o',actionlistdecreasek(:,1),actionlistdecreasek(:,2),'o','MarkerSize',4)
hold on
plot(kappamesh, 0*kappamesh -2*betaj+pi^2/2, kappamesh,-2*betaj*kappamesh/sqrt(lz)+pi^2/2/lz,'LineWidth',1)
     set(gca,'FontSize',16)

titlestr= append("$\beta \mathbf{J}$ = " , num2str(betaj, '%.4g'), ", $l_{nz} =$ " , num2str(l, '%.4g'),", $l_{zz}$ = " , num2str(lz, '%.4g')) 
 title(titlestr,'Interpreter','latex',FontSize=20)
  xlabel('$\kappa$','Interpreter','latex',FontSize=20)
  ylabel('Action', 'Interpreter','latex', FontSize=20)
  lgstr1 ='increasing $\kappa$ numerical';
  lgstr2='decreasing $\kappa$ numerical';
  lgstr3 = 'increasing $\kappa$ analytic approx.';
  lgstr4= 'decreasing $\kappa$ analytic approx.';
  legend(lgstr1,lgstr2,lgstr3,lgstr4,'Interpreter','latex','Location', 'southwest',FontSize=18)
grid on
 %actionlistdecreasek1000 =  actionlistdecreasek;
 %actionlistincreasek1000=actionlistincreasek;
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

