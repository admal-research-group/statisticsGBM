% This code runs grain growth simulation with simple anisotropic grain boundary
% characters. In particular, grains are grouped into 3 types, A,B,C. Grain boundaries
% between different groups have larger grain boundary energy and higher
% mobilities, while the same group has lower energy and mobiltiy.


clear
close all 

addpath(genpath(pwd)) % add all subfolders to the path

%% Hyper parameters and Global variables 

FinalTime = 50; % Final time-step to be simulated 
N=1000; % The number of grains
dims=[1000 1000]; 
dt = 5/(dims(1)^2); % suggested time step size 

Rgrains = 10; % number of pixel to outgrow, size of buffer region for each grain 
Rfamilies = 30; % minimum pixel distance to have two grains in sthe ame family 
I = sqrt(-1); % Imaginary number i.
w=exp(I*2*pi/dims(1)); % nth root of unity.x
[x,y] = meshgrid(1:dims(1),1:dims(1)); x = x'; y = y';

%% WORK SPACES 
Z = -ones(dims(1),dims(2)); % workspace needed for loc_levset_to_volfluid.c.
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

% Create Grain 
[grains,cluster] = EModule_initialdoubleVoronoidata2d(N,dims,Rgrains,rn,FactorNc);
% Grain IDs
ID = (1:1:size(grains,1))';

% Create Orientations 
ori=zeros(N,1); 

% First, we generate integer type ori, which will be used to preprare Gaussian kernels 
for i=1:length(ori)
    
 if(mod(i,3)==0 )
 ori(i) = randi([0,3]); % from 0 to 5 to be red (Group A) 
 elseif(mod(i,3)==1)
 ori(i) = 12 +  3 * randi([0,3]); %from 15 to 20, to be green (Group B)
 elseif(mod(i,3)==2)
 ori(i) = 62 +  3 * randi([0,3]); % from 25 to 30, to be red (Group C) 
 end
      
end

ori = ori * pi/180 ; % make radian 

% Determine the width of two Gaussian Kernels 

distinctOri= unique(ori); 
kernelOption = 2;
[alpha, beta] = EModule_kernel_widths2d_ABC(distinctOri,kernelOption);
 
% Create Gaussian Kernels 
n=dims(1); 
KERNELalpha = exp(-dt*alpha*n*n*(4-w.^(x-1)-w.^(1-x)-w.^(y-1)-w.^(1-y)));
KERNELbeta = exp(-dt*beta*n*n*(4-w.^(x-1)-w.^(1-x)-w.^(y-1)-w.^(1-y)));

disp('alpha'); disp(alpha);
disp('beta'); disp(beta); 


% After kerenls are prepared, we regenerate ori in float type
for i=1:length(ori)
    
 if(mod(i,3)==0 )
 ori(i) = 3 * rand(); % from 0 to 5 to be red
 elseif(mod(i,3)==1)
 ori(i) = 13 +  3 * rand(); %from 15 to 20, to be green
 elseif(mod(i,3)==2)
 ori(i) = 62 +  3 * rand(); % from 25 to 30, to be red
 end

end

ori = ori * pi/180 ; % make radian

% Time marching
ori_in_deg = ori * 180/pi; 

fileName1 = './ABCToyMobil_';
for t = 1:FinalTime % Main time iteration starts.
    
    status = append(int2str(t),' Time step is running');
    disp(status); 
    
    %EModule_savegrainsfigure_withLabelRGB_ori(grains,dims,ori,fileName1,t); 

	[families,famgrains] = GModule_grains2families2d(grains,Rfamilies,dims,WORKSPACE,work_x, work_y, work_bx, work_by, work_candidate_x, work_candidate_y); 
    Nf = size(families,1); % Number of families.
   
    %if( t>10 && mod(t,50)==0)
  	%filename='./ABCToy_'; 
	%filename= append(filename,int2str(t)); 
    %save(filename,'grains','dims','-v7.3')
    %end

     
    % Convolution     
    for k=1:Nf % Loop over families.
        
        % convolute on each family 
        cval = GModule_convolvefamily2d(families{k,1},families{k,2},dims,KERNELalpha, KERNELbeta, Z);
        numberofgrains = size(famgrains{k},1); % the number of grains in this family.
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

    option=2; angBrandon =30; % option 2 : step mobility 
    CEModule_updatelevelsetdata2d_ABC(presence,grains,ID,ori_in_deg,alpha,beta,angBrandon,option);
    
    % Thresholding:
    grains = GModule_Thresholding(grains,Rgrains,dims,Z,WORKSPACE,...
            work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y);
    
end

