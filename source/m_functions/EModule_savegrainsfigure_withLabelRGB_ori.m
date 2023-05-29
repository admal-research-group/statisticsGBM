function EModule_savegrainsfigure_withLabelRGB_ori(grains,dims,ori,fileName, fileNumber)
% outputs grains configuration into an *png* image file as a part of  a time-series figure. While the color scales of grains are determined by orientation values, the specific rule must be rewritten as needed. 
%
% :param grains: data structure grains 
% :param dim: the size of the system
% :param ori: orientations of grains in 1d array, which determines the color scales.  
% :param fileName: a string for output file name 
% :param fileNumber: the number of figure in time-series
% :return: an image file. e.g., fileName_000%d.png 

Ng = size(grains,1); % Number of grains.

% matlab rgb colors takes values in between [0,1]
% the following code will well resolve large polycrystals....

u1 = zeros(dims);
u2 = zeros(dims);
u3 = zeros(dims); 


for k=1:Ng % Loop over grains.
  ind = grains{k,1}; % Pixels in a nhd. of the grain.
  val = grains{k,2}; % Level set values. 
  ind2 = ind(val>0); % Pixels in the interior of grain.

if( ori(k) >0.0 && ori(k)*180/pi <=5) 
       colorscale = (ori(k)*180/pi) / (3*2);
       u1(ind2) =  0.5 + colorscale; %red
elseif( ori(k)*180/pi <=20)        
       colorscale = (ori(k)*180/pi-13) / (3*2) ;
       u2(ind2) = 0.5 + colorscale;  %green
elseif( ori(k)*180/pi <=30)        
       colorscale = (ori(k)*180/pi-26) / (3*2) ;
       u3(ind2) = 0.5 + colorscale; %emerald
       u2(ind2) = 0.5 + colorscale;
elseif( ori(k)*180/pi <=50) 
       colorscale = (ori(k)*180/pi-39) / (3*2) ; 
       u1(ind2) = 0.5 + colorscale;% yellow
       u2(ind2) = 0.5 + colorscale;% 
elseif( ori(k)*180/pi <=70) 
       colorscale = (ori(k)*180/pi-63) / (3*2) ;  
       u3(ind2) = 0.5 + colorscale; % blue
end
 
end


data(:,:,1)=u1; 
data(:,:,2)=u2; 
data(:,:,3)=u3; 

fig=figure; 


image(data);
axis square;
axis off;

% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];

formatSpec='%04d'; 
formatSpec2='.png'; 

ax = gcf;
FnameString = append(fileName,num2str(fileNumber, formatSpec)) ; 
FnameString = append(FnameString,formatSpec2) ; 

%print(fig,FnameString,'-dpng','-r0')

exportgraphics(ax,FnameString,'Resolution',300) 

close(fig); 


end
