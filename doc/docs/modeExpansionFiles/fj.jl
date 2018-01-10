using HDF5;
using PyPlot;

@printf("A\n");
file=h5open("fj1_L0_fluxA.h5");
exA=read(file,"ex_0.r") + im*read(file,"ex_0.i");
eyA=read(file,"ey_0.r") + im*read(file,"ey_0.i");
hxA=read(file,"hx_0.r") + im*read(file,"hx_0.i");
hyA=read(file,"hy_0.r") + im*read(file,"hy_0.i");
close(file)

@printf("B\n");
file=h5open("fj1_L0_fluxB.h5");
exB=read(file,"ex_0.r") + im*read(file,"ex_0.i");
eyB=read(file,"ey_0.r") + im*read(file,"ey_0.i");
hxB=read(file,"hx_0.r") + im*read(file,"hx_0.i");
hyB=read(file,"hy_0.r") + im*read(file,"hy_0.i");
close(file)

@printf("1\n");
file=h5open("fj1_L0_mode1.h5");
ex1=read(file,"ex.r") + im*read(file,"ex.i");
ey1=read(file,"ey.r") + im*read(file,"ey.i");
hx1=read(file,"hx.r") + im*read(file,"hx.i");
hy1=read(file,"hy.r") + im*read(file,"hy.i");
close(file)

@printf("2\n");
file=h5open("fj1_L0_mode2.h5");
ex2=read(file,"ex.r") + im*read(file,"ex.i");
ey2=read(file,"ey.r") + im*read(file,"ey.i");
hx2=read(file,"hx.r") + im*read(file,"hx.i");
hy2=read(file,"hy.r") + im*read(file,"hy.i");
close(file)

#=
@printf("3\n");
file=h5open("fj1_L0_mode3.h5");
ex3=read(file,"ex.r") + im*read(file,"ex.i");
ey3=read(file,"ey.r") + im*read(file,"ey.i");
hx3=read(file,"hx.r") + im*read(file,"hx.i");
hy3=read(file,"hy.r") + im*read(file,"hy.i");
close(file)

@printf("4\n");
file=h5open("fj1_L0_mode4.h5");
ex4=read(file,"ex.r") + im*read(file,"ex.i");
ey4=read(file,"ey.r") + im*read(file,"ey.i");
hx4=read(file,"hx.r") + im*read(file,"hx.i");
hy4=read(file,"hy.r") + im*read(file,"hy.i");
close(file)
=#

@printf("C\n");
file=h5open("fj1_L0_fluxC.h5");
eyC=read(file,"ey_0.r") + im*read(file,"ey_0.i");
ezC=read(file,"ez_0.r") + im*read(file,"ez_0.i");
hyC=read(file,"hy_0.r") + im*read(file,"hy_0.i");
hzC=read(file,"hz_0.r") + im*read(file,"hz_0.i");
close(file)

nzA=2;
nzB=1;
nz1=1;
nz2=1;

LX=3.0;
LY=3.0;
clf();
subplot(4,4,1)
imshow(abs(exA[nzA,:,:])); colorbar()
subplot(4,4,2)
imshow(abs(eyA[nzA,:,:])); colorbar()
subplot(4,4,3)
imshow(abs(hxA[nzA,:,:])); colorbar() 
subplot(4,4,4)
imshow(abs(hyA[nzA,:,:])); colorbar() 

subplot(4,4,5)
imshow(abs(exB[nzB,:,:])); colorbar()
subplot(4,4,6)
imshow(abs(eyB[nzB,:,:])); colorbar()
subplot(4,4,7)
imshow(abs(hxB[nzB,:,:])); colorbar()
subplot(4,4,8)
imshow(abs(hyB[nzB,:,:])); colorbar()

subplot(4,4,9)
imshow(abs(ex1[nz1,:,:])); colorbar()
subplot(4,4,10)
imshow(abs(ey1[nz1,:,:])); colorbar()
subplot(4,4,11)
imshow(abs(hx1[nz1,:,:])); colorbar()
subplot(4,4,12)
imshow(abs(hy1[nz1,:,:])); colorbar()

subplot(4,4,13)
imshow(abs(ex2[nz2,:,:])); colorbar()
subplot(4,4,14)
imshow(abs(ey2[nz2,:,:])); colorbar()
subplot(4,4,15)
imshow(abs(hx2[nz2,:,:])); colorbar()
subplot(4,4,16)
imshow(abs(hy2[nz2,:,:])); colorbar()
