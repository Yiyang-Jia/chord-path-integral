%lz=0;
%l=1;



clear all
global  betaj k v l lz

l=1;
lz=1.5;
k=0.72; %at fixed k, large betaj tend to bring you to p-spin phase.
% 2nd order critical value (betja =13.5, k =0.457)
betajlarge = 1000;
betajsmall = 60;
%samplebetaj=30;%look at the solution(s) at this point
betaj = betajsmall; %start from SYK phase and increase betaj

syms f(x)
f(x) = pi*x/cos(pi* x/2) - betaj;
v = eval(vpasolve(f==0, x, [0 1]));

xmesh=linspace(-0.5, 0.5, 553);
 options = bvpset('RelTol',1e-4, 'Nmax', 5000);
%solinitSmallbeta = bvpinit(xmesh, [-0.5; 0;-3;0]);
solinitSmallbeta = bvpinit(xmesh, @exactk0);
solSmallbeta= bvp4c(@bvpfun, @bcfun, solinitSmallbeta, options);

actionlistincreasebetaj = [betaj, action(solSmallbeta)];  %values of action, with k increasing from 0 
soliterateInc = solSmallbeta;

spacing = 2;
 npoints = floor((betajlarge - betajsmall)* (1/spacing));

for n = 1: npoints
         lastwarn('') % Clear last warning message
          warnIncreasebetaj="";
    errIncreasebetaj = "";
      betaj = betaj + spacing; 

  try
    solSmallbeta =bvpinit(soliterateInc.x, @(u) interp1(soliterateInc.x, transpose(soliterateInc.y), u));%https://www.mathworks.com/matlabcentral/answers/142088-initial-guess-for-bvp4c
    soliterateInc= bvp4c(@bvpfun, @bcfun, solinitSmallbeta, options);
    actionlistincreasebetaj = [actionlistincreasebetaj; betaj action(soliterateInc)]
    % if abs(betaj - samplebetaj)<0.3
    %     solSampleInc = soliterateInc;
    % end
    [warnMsg, warnId] = lastwarn;
   
     if ~isempty(warnMsg)%break the loop is there's a warning message (poor convergence when solving ode)
     warnIncreasebetaj = "Warning: increasing betaj sequence, betaj= "+ betaj + ". " + warnMsg;
       break 
     end
  catch err
       errIncreasebetaj = "Error: increasing betaj sequence, betaj= "+ betaj + ". " + err.message;
      break
  end 
end
actionlistincreasebetaj(end,:) = []; %remove the last result when the loop breaks, this is a poorly converged result

%% 

betaj = betajlarge;

xmesh=linspace(-0.5, 0.5, 553);
solinitLargebeta = bvpinit(xmesh, @exactk1lz);
%solinitLargebeta = bvpinit(xmesh, @exactk1);
%solinitLargebeta = bvpinit(xmesh, [-3;0;-1;0]);
solLargebeta = bvp4c(@bvpfun, @bcfun, solinitLargebeta,options);
actionlistdecreasebetaj = [betaj, action(solLargebeta)];  %actions, with k decreasing from 1

soliterateDec =solLargebeta;


for n = 1: npoints
  
       lastwarn('') % Clear last warning message
       warnDecreasebetaj = "";
       errDecreasebetaj = "" ;
      betaj = betaj - spacing; 
      
      try
      solLargebeta =bvpinit(soliterateDec.x, @(u) interp1(soliterateDec.x, transpose(soliterateDec.y), u));%https://www.mathworks.com/matlabcentral/answers/142088-initial-guess-for-bvp4c
      soliterateDec= bvp4c(@bvpfun, @bcfun, solinitLargebeta,options);
      actionlistdecreasebetaj = [actionlistdecreasebetaj; betaj action(soliterateDec)]
        % if abs(betaj - samplebetaj)<0.3
        %   solSampleDec = soliterateDec;
        % end
       [warnMsg, warnId] = lastwarn;
       
      if ~isempty(warnMsg)%break the loop is there's a warning message (poor convergence when solving ode)
     warnDecreasebetaj = "Warning: decreasing betaj sequence, betaj= "+ betaj + ". " + warnMsg;
       break 
      end
      
  catch err
      errDecreasebetaj = "Error: decreasing betaj sequence, betaj= "+ betaj + ". " + err.message;
      break
   end
      if betaj <0
          break
      end
end
 actionlistdecreasebetaj(end,:) = [];
 
 fprintf( '%s\n',warnIncreasebetaj)
   fprintf('%s\n', errIncreasebetaj)

 fprintf('%s\n',warnDecreasebetaj)
 fprintf('%s\n', errDecreasebetaj)
% 

%% 
figure;
 plot(actionlistincreasebetaj(:,1),actionlistincreasebetaj(:,2),actionlistdecreasebetaj(:,1),actionlistdecreasebetaj(:,2))
  title("Action as a function of \beta*J  (k = " + num2str(k, '%.2f')+", l = " + num2str(l, '%.2f')+", l_z = " + num2str(lz, '%.2f')+")" )
  xlabel('\beta*J ')
  ylabel('S')

  %% 
   figure;
  %  subplot(1,2,1)
  % plot(solSampleInc.x, solSampleInc.y(1,:),solSampleInc.x, solSampleInc.y (3,:))
  %  subplot(1,2,2)
  %  plot(solSampleDec.x, solSampleDec.y(1,:),solSampleDec.x, solSampleDec.y (3,:))
subplot(1,2,1)
  plot(soliterateInc.x, soliterateInc.y(1,:),soliterateInc.x, soliterateInc.y (3,:))
  subplot(1,2,2)
  plot(soliterateDec.x, soliterateDec.y(1,:),soliterateDec.x, soliterateDec.y (3,:))

   %% 
   
   % transitionPts = findcritk(actionlistincreasebetaj,actionlistdecreasebetaj)
 % dlmwrite('tranpts.txt',transitionPts,'-append');
 %in the format of [1/betaj  kc]



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
