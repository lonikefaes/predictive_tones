function Dr = upsampling2_5(D,interp_factor)

[nx0,ny0,nz0] = size(D);

[X0 Y0 Z0]    = meshgrid(0:ny0-1,0:nx0-1,0:nz0-1);
[Xr Yr Zr]    = meshgrid(linspace(0,ny0-1,ny0*interp_factor/2),...
                         linspace(0,nx0-1,nx0*interp_factor/2),...
                         linspace(0,nz0-1,nz0*interp_factor/2));

Dr = interp3(X0,Y0,Z0,D,Xr,Yr,Zr,'nearest');