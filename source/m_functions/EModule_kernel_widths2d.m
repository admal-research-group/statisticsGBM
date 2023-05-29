function [alpha,beta] = EModule_kernel_widths2d(distinctOri,energy_curve,option)

% Estimates the characteristic widths of two Gaussian kernels using the external GB energy data
%
% :param distinctOri: an 1d array of distinct orientations (non-repeated)
% :param energy_curve: an 1d array of distinct orientations (non-repeated)
% :param option: option for mobility, option 1: unit moblity, option 2: high mobility for high misorientation angle grain boundaries
% :returns: ``[alpha, beta]`` characteristic widths of two Gaussian kernels

N = length(distinctOri); % returns the number of grain

if(option > 0)    
  S = ones(N,N)-eye(N,N); % Surface tension matrix 
  
    for i=1:N
     for j=i+1:N
      % convert to degree...   
      oriA=distinctOri(i)* 180/pi; oriB=distinctOri(j)*180/pi;
      
      misorient= abs(oriA-oriB);
      energy = interp1(energy_curve(:,1),energy_curve(:,2),misorient);
      S(i,j) = energy; 
      S(j,i) = S(i,j); 
     end
    end
  
  % High mobility for high angle grain boundaries 
  % Low mobility for low angle grain boundaries 
  if(option==2)
    iM = zeros(N,N) ; % inverseMobility     
    % Mobility function parameters 
    thetaMax = 10; 
    
    for i=1:N
     for j=i+1:N
      % convert to degree...   
      oriA=distinctOri(i)* 180/pi; oriB=distinctOri(j)*180/pi;
      misorient= abs(oriA-oriB); 
       if(misorient < thetaMax)
          mobility = 0.905;
       else
          mobility = 1.0; 
       end
       
      iM(i,j) = 1.0/mobility;
      iM(j,i) = iM(i,j); 
      
     end
    end
    
  end   
  
    
    
elseif(option==-1)

  minang = misorientation_angle2d(distinctOri);
  S = ones(N,N)-eye(N,N); % Surface tension matrix 
  angBrandonRad = 30 * pi/180;  
  select = minang<=angBrandonRad;
  S(select) = minang(select)/angBrandonRad.*(1-log(minang(select)/angBrandonRad));
  S(1:N+1:N^2) = 0;
    
end


% Read-Shockley with Bradon angle %
% sufrace_tension2d function looks as follow %

J = eye(N)-1/N*ones(N,1)*ones(1,N);
eigS = eig(J*S*J);
[~,i] = min(abs(eigS));
eigS(i) = [];

absoluteEigenS=-abs(eigS);  
mineigS = min(absoluteEigenS);
maxeigS = max(absoluteEigenS);


% Get max and min of reciprocal M 
if (option==1 || option==-1) % All mobilities = 1
 maxeigreciprocalM = -1;
 mineigreciprocalM = -1;
elseif(option ==2) 
 eigiM = eig(J*iM*J);
 [~,i] = min(abs(eigiM)); % This get rid of unstable value 
 eigiM(i) = [];
 
 assert( max(eigiM) < 0 ) ; % for numberical stability 
 absoluteEigeniM=-abs(eigiM);  % These are all negative 
 maxeigreciprocalM = max(absoluteEigeniM); 
 mineigreciprocalM = min(absoluteEigeniM); 
end

alpha = mineigS/maxeigreciprocalM;
beta = maxeigS/mineigreciprocalM;
 
M = ones(N)-eye(N);
aux = S.*M;
aux(1:N+1:N^2) = -Inf;
alpha = max(alpha,max(aux(:)));
aux(1:N+1:N^2) = Inf;
beta = min(beta,min(aux(:)));

if sum(eigS > 0)
   error('The matrix \sigma is not conditionally negative semidefinite.')
end

end

