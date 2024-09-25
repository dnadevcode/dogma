
barLong = zeros(1,300)
barLong(150) = 1;
PSF_WIDTH_PIXELS  = 300/110;
barLong = imgaussfilt(barLong,PSF_WIDTH_PIXELS,'Padding','circular');

% [find(barLong>0,1,'Last')-find(barLong>0,1,'First')]*max(barLong)/

max(barLong)/sum(barLong)
(1/sqrt(2*pi)*1/PSF_WIDTH_PIXELS)/2