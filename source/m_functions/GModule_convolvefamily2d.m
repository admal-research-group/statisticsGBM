function cval = GModule_convolvefamily2d(ind,val,dims,KERNELalpha, KERNELbeta, Z)
% Convolves the characteristic function of a family (a union of grains)
% whose level set representation is stored in *ind* and *val* input vars using the kernel *KERNEL*.
%
% :param ind: grid indices that belongs to the family 
% :param val: the level-set value of the characteristic function of the family 
% :param dim: the size of the system
% :param Kernelalpha: Gaussian Kernel with characteristic width alpha, size of *dim*
% :param Kernelbeta: Gaussian Kernel with characteristic width beta, size of *dim*
% :param Z: Work space, size of *dim*
% :return cval: Convolution values at those grid indices

% Convert level set data to volume of fluid representation:
[x,y] = ind2sub(dims,ind);
vf = CGModule_loc_levset_to_volfluid(int32(x),int32(y),val,Z);

% Carry out the convolution:
Kalphau = zeros(dims);
Kalphau(ind) = vf; % Characteristic function of the union of grains.
Kalphau = real(ifft2(fft2(Kalphau).*KERNELalpha));
Kbetau = zeros(dims);
Kbetau(ind) = vf; % Characteristic function of the union of grains.
Kbetau = real(ifft2(fft2(Kbetau).*KERNELbeta));
cval{1} = Kalphau; % Convolution values, returned.
cval{2} = Kbetau;

end

