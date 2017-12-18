using HDF5;
using PyPlot;

file=h5open("dj_fluxA.h5");
exA=read(file,"ex_0.r") + im*read(file,"ex_0.i");
eyA=read(file,"ey_0.r") + im*read(file,"ey_0.i");
hxA=read(file,"hx_0.r") + im*read(file,"hx_0.i");
hyA=read(file,"hy_0.r") + im*read(file,"hy_0.i");
close(file)

file=h5open("dj_fluxB.h5");
exB=read(file,"ex_0.r") + im*read(file,"ex_0.i");
eyB=read(file,"ey_0.r") + im*read(file,"ey_0.i");
hxB=read(file,"hx_0.r") + im*read(file,"hx_0.i");
hyB=read(file,"hy_0.r") + im*read(file,"hy_0.i");
close(file)

file=h5open("dj_mode1.h5");
ex1=read(file,"ex.r") + im*read(file,"ex.i");
ey1=read(file,"ey.r") + im*read(file,"ey.i");
hx1=read(file,"hx.r") + im*read(file,"hx.i");
hy1=read(file,"hy.r") + im*read(file,"hy.i");
close(file)

file=h5open("dj_mode2.h5");
ex2=read(file,"ex.r") + im*read(file,"ex.i");
ey2=read(file,"ey.r") + im*read(file,"ey.i");
hx2=read(file,"hx.r") + im*read(file,"hx.i");
hy2=read(file,"hy.r") + im*read(file,"hy.i");
close(file)

file=h5open("dj_mode3.h5");
ex3=read(file,"ex.r") + im*read(file,"ex.i");
ey3=read(file,"ey.r") + im*read(file,"ey.i");
hx3=read(file,"hx.r") + im*read(file,"hx.i");
hy3=read(file,"hy.r") + im*read(file,"hy.i");
close(file)

file=h5open("dj_mode4.h5");
ex4=read(file,"ex.r") + im*read(file,"ex.i");
ey4=read(file,"ey.r") + im*read(file,"ey.i");
hx4=read(file,"hx.r") + im*read(file,"hx.i");
hy4=read(file,"hy.r") + im*read(file,"hy.i");
close(file)

file=h5open("dj_fluxC.h5");
exC=read(file,"ex_0.r") + im*read(file,"ex_0.i");
ezC=read(file,"ez_0.r") + im*read(file,"ez_0.i");
hxC=read(file,"hx_0.r") + im*read(file,"hx_0.i");
hzC=read(file,"hz_0.r") + im*read(file,"hz_0.i");
close(file)

sA=zeros(size(exB));
sB=zeros(size(exB));
s1=zeros(size(ex1));
s2=zeros(size(ex2));
s3=zeros(size(ex3));
s4=zeros(size(ex4));
for m=1:size(sA,1), n=1:size(sA,2)
  sA[m,n] = real( conj(exA[m,n])*hyA[m,n] - conj(eyA[m,n])*hxA[m,n] );
  sB[m,n] = real( conj(exB[m,n])*hyB[m,n] - conj(eyB[m,n])*hxB[m,n] );
  s1[m,n] = real( conj(ex1[m,n])*hy1[m,n] - conj(ey1[m,n])*hx1[m,n] );
  s2[m,n] = real( conj(ex2[m,n])*hy2[m,n] - conj(ey2[m,n])*hx2[m,n] );
  s3[m,n] = real( conj(ex3[m,n])*hy3[m,n] - conj(ey3[m,n])*hx3[m,n] );
  s4[m,n] = real( conj(ex4[m,n])*hy4[m,n] - conj(ey4[m,n])*hx4[m,n] );
end


clf();
subplot(3,4,1)
imshow(abs(exA)); title("ex(A)"); xlabel("y"); ylabel("|ex|");
colorbar()

subplot(3,4,2)
imshow(abs(eyA)); title("ey(A)"); xlabel("y"); ylabel("|ey|");
colorbar()

subplot(3,4,3)
imshow(abs(hxA)); title("hx(A)"); xlabel("y"); ylabel("|hx|");
colorbar()

subplot(3,4,4)
imshow(abs(hyA)); title("hy(A)"); xlabel("y"); ylabel("|hy|");
colorbar()

subplot(3,4,5)
imshow(abs(exB)); title("ex(B)"); xlabel("y"); ylabel("|ex|");
colorbar()

subplot(3,4,6)
imshow(abs(eyB)); title("ey(B)"); xlabel("y"); ylabel("|ey|");
colorbar()

subplot(3,4,7)
imshow(abs(hxB)); title("hx(B)"); xlabel("y"); ylabel("|hx|");
colorbar()

subplot(3,4,8)
imshow(abs(hyB)); title("hy(B)"); xlabel("y"); ylabel("|hy|");
colorbar()

subplot(3,4,9)
imshow(abs(ex1)); title("ex(1)"); xlabel("y"); ylabel("|ex|");
colorbar()

subplot(3,4,10)
imshow(abs(ey1)); title("ey(1)"); xlabel("y"); ylabel("|ey|");
colorbar()

subplot(3,4,11)
imshow(abs(hx1)); title("hx(1)"); xlabel("y"); ylabel("|hx|");
colorbar()

subplot(3,4,12)
imshow(abs(hy1)); title("hy(1)"); xlabel("y"); ylabel("|hy|");
colorbar()
