function newgrains = GModule_Thresholding(grains,Rgrains,dims,Z,WORKSPACE,work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y)
% Obtain a new grain using thresholding of the comparion functions
%
% :param grains: data structure grains
% :param Z: pre-allocated workspace 
% :param WORKSPACE: pre-allocated workspace
% :param work_x: pre-allocated workspace, same for others 
% :returns: new grains some outputs

newgrains = grains; % prepare data space 

Ng = size(grains,1);

for k=1:Ng % Loop over grains.
  ind = grains{k,1}; % Pixels within a nhd. of grain.
  val = grains{k,2}; % Lev. set. vals. at those pixels.
  cval1 = grains{k,3}; % Convolution vals. at those pixels.
  cval3 = grains{k,4}; % Convolution vals. at those pixels.
        
  Z(ind) = val;      % Lev. set. representation on grid.
  posind = ind(val>0); % Pixels in the interior of grain.
  [x,y] = ind2sub(dims,posind);
  [x2,y2] = CGModule_pgrow3(int32(x),int32(y),Rgrains,WORKSPACE,...
            work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y); % Dilation.
  ind2 = sub2ind(dims,x2,y2);
  val2 = Z(ind2); % Level set vals.
  Z(ind2) = -1; % Reset Z back to all -1's.
  Z(ind) = cval1 - 1; % Convolution values - 1.
  cval2 = Z(ind2); % Convolution vals - 1.
  Z(ind2) = -1;
  Z(ind) = cval3 - 1; % Convolution values - 1.
  cval4 = Z(ind2); % Convolution vals - 1.
  Z(ind2) = -1;
  
  newgrains{k,1} = ind2;   % Refresh grain's data structure.
  newgrains{k,2} = val2;   % Ditto.
  newgrains{k,3} = cval2 + 1; % Ditto.
  newgrains{k,4} = cval4 + 1; % Ditto.
end % (for k). Loop over grains ends.

end