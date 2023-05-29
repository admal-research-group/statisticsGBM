% This code runs grain growth simulation with [110] symmetric tilt grain
% boundary energies. Orientation values of grains are sampled from
% sub-domain region with misorientation angle [0,70.6] degree from a reference state 
% Also, mobilities for small angle grain boundaries are set to lower than that of
% high angle grain boundaries. 

clear
close all 

addpath(genpath(pwd)) % add all subfolders to the path

%% Hyper parameters 
FinalTime = 20;

N=1000;
dims=[1000 1000]; 
dt = 20/(dims(1)^2); % suggested time step size 

Rgrains = 10; % number of pixel to outgrow 
Rfamilies = 30; % minimum pixel distance to have two grains in same family 
I = sqrt(-1); % Imaginary number i.
w=exp(I*2*pi/dims(1)); % nth root of unity.x
[x,y] = meshgrid(1:dims(1),1:dims(1)); x = x'; y = y';

% Read the external data file
CVE = importdata("data/FCC_110STGB.txt"); % external GB energy data file 

%% WORK SPACES 
Z = -ones(dims(1),dims(2)); % workspace needed in loc_levset_to_volfluid.c.
WORKSPACE = int32(-ones(dims)); % workspace needed for pgrow
work_x = int32(zeros(prod(dims),1));
work_y = int32(zeros(prod(dims),1));
work_bx = int32(zeros(prod(dims),1));
work_by = int32(zeros(prod(dims),1));
work_candidate_x = int32(zeros(prod(dims),1));
work_candidate_y = int32(zeros(prod(dims),1));

%% Initialization
rn=1; % random seed number  
rng(rn); % Control random number generator
Nc=30; 
FactorNc = fix(N/Nc);

[grains, cluster] = EModule_initialdoubleVoronoidata2d(N,dims,Rgrains,rn,FactorNc);

ori= zeros(N,1); 

for i=1:length(ori)
    
 if(mod(i,5)==0 )
 ori(i) = 3 * rand(); % from 0 to 5
 elseif(mod(i,5)==1)
 ori(i) = 13 +  3 * rand(); %
 elseif(mod(i,5)==2)
 ori(i) = 26 +  3 * rand(); % 
 elseif(mod(i,5)==3)
 ori(i) = 39 +  3 * rand(); % 
 elseif(mod(i,5)==4)
 ori(i) = 63 +  3 * rand(); % 
 end
      
end

ori_in_deg = ori; 
ori = ori * pi/180; % make radian 

% Now, we need to prepare the Gaussian kernels
% prepare temporary list of ori values for kernel estimation
ori_for_kernel = zeros(N,1); 
ID = (1:1:size(grains,1))';

for i=1:length(ori_for_kernel) 

 if(mod(i,5)==0 )
 ori_for_kernel(i) = 0 +   0.5*randi([0,6]); % from 0 to 3
 elseif(mod(i,5)==1)
 ori_for_kernel(i) = 13 +  0.5*randi([0,6]); %from 13 to 16
 elseif(mod(i,5)==2)
 ori_for_kernel(i) = 26 +  0.5*randi([0,6]); % from 26 to 29
 elseif(mod(i,5)==3)
 ori_for_kernel(i) = 39 +  0.5*randi([0,6]); % from 39 to 42 
 elseif(mod(i,5)==4)
 ori_for_kernel(i) = 63 +  0.5*randi([0,6]); % from 63 to 66
 end
      
end

ori_for_kernel = ori_for_kernel * pi/180 ; % make radian 
distinctOri= unique(ori_for_kernel); 

% Build Kernel for unit mobility 
option=2; % Step Mobility 
[alpha, beta] = EModule_kernel_widths2d(distinctOri,CVE,option);

% Create Gaussian Kernels 
n=dims(1); 
KERNELalpha = exp(-dt*alpha*n*n*(4-w.^(x-1)-w.^(1-x)-w.^(y-1)-w.^(1-y)));
KERNELbeta = exp(-dt*beta*n*n*(4-w.^(x-1)-w.^(1-x)-w.^(y-1)-w.^(1-y)));

disp('alpha'); disp(alpha);
disp('beta'); disp(beta); 


fileName1 = './HAGBpl_';

for t = 1:FinalTime % Main time iteration starts.
    
    status = append(int2str(t),' Time step is running');
    disp(status); 
    
    Module_DataProcess_savegrainsfigure_withLabelRGB_ori(grains,dims,ID,ori,fileName1,t); 
    
	[families,famgrains] = GModule_grains2families2d(grains,Rfamilies,dims,WORKSPACE,work_x, work_y, work_bx, work_by, work_candidate_x, work_candidate_y);
    Nf = size(families,1); % Number of families.
   
    % save current grain configuration in mat data file 
    %if( mod(t,500)==0)
    %filename='./HAGB_pl_'; 
    %filename= append(filename,int2str(t)); 
    %save(filename,'grains','dims','ori','-v7.3')
    %end

    % Convolving families     
    for k=1:Nf % Loop over families.
        
        % convolute on each family 
        cval = GModule_convolvefamily2d(families{k,1},families{k,2},dims,KERNELalpha, KERNELbeta, Z);
        numberofgrains = size(famgrains{k},1); % Number of grains in this family.
        listofgrains = famgrains{k};           % Column vector of grain indices.
        
        % redistribute convolution values to the grains contained in this family:
        for ell = 1:numberofgrains % Loop over grains contained in this family.
            label = listofgrains(ell);
            ind = grains{label,1};
            cvalaux = cval{1}; 
            grains{label,3} = cvalaux(ind); % Read off and record the convolution vals.
            cvalaux = cval{2};
            grains{label,4} = cvalaux(ind); 
        end % (for ell) Loop over grains  ends.
    end % (for k) Loop over families ends.
    
    % REDISTRIBUTION STEP:
    presence = CGModule_get_nhd_grains2d(grains,dims(1)*dims(2));

    option=2; angBrandon =30; 
    CEModule_updatelevelsetdata2d_CVE(presence,grains,ID,ori_in_deg,alpha,beta,angBrandon,option,CVE);
   
    % Thresholding:
    grains = GModule_Thresholding(grains,Rgrains,dims,Z,WORKSPACE,...
            work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y);
end

