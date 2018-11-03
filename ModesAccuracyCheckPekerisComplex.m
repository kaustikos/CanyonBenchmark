function  err = ModesAccuracyCheckPekerisComplex(wnum,MediaParams,freq)


if isstruct(MediaParams)
    MP = MediaParams.LayersData;
else
    MP = MediaParams;
end;



omeg = 2*pi*freq;
H = MP(2,1);
cw = MP(2,2);
cb = MP(2,3);
rhow = MP(2,4);
rhob = MP(2,5);


eta = 1/(40*pi*log10(exp(1)));    
betab = MP(2,7);
cb = cb/(1+1i*eta*betab);


kvw = sqrt( (omeg/cw)^2 - wnum.^2 );
kvb = sqrt(  wnum.^2 - (omeg/cb)^2 );


err = tan( kvw*H ) + rhob*kvw./( rhow*kvb );