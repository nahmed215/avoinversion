function V=vel_smoother(c, pad, perc1, perc2, tap)
[Nx2, Nz2]=size(c);
if(nargin <5)
    tap=1;
end
cleft=c(:,1)*ones(1, pad)*tap;
cright=c(:,end)*ones(1, pad)*tap;
c=[cleft, c, cright];
ctop=ones(pad,1)*c(1, :)*tap;
cbottom=ones(pad,1)*c(end,:)*tap;
c=[ctop;c;cbottom];
clear ctop, clear cbottom, clear cleft, clear cright;

[Nx, Nz]=size(c);
N1=2^nextpow2(Nx);
N2=2^nextpow2(Nz);

L1=floor(perc1*N1);
L2=ceil((N1-L1)/2);
L1=N1-L2*2;
H1=fftshift([zeros(1, L2) 0.5*(1-cos(2*pi*[0:((L1-1))]/(L1-1))) zeros(1, L2)]);
 
L1=floor(perc2*N2);
L2=ceil((N2-L1)/2);
L1=N2-L2*2;
H2=fftshift([zeros(1, L2) 0.5*(1-cos(2*pi*[0:((L1-1))]/(L1-1))) zeros(1, L2)]);

H=H1.'*H2;
c2=real(ifft2(fft2(c, N1, N2).*H));

V=c2(1+pad:Nx2+pad, 1+pad:Nz2+pad);



