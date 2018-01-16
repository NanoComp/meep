using HDF5;
using PyPlot;

file=h5open("pj_fluxA.h5");
exA=read(file,"ex_0.r") + im*read(file,"ex_0.i");
eyA=read(file,"ey_0.r") + im*read(file,"ey_0.i");
hxA=read(file,"hx_0.r") + im*read(file,"hx_0.i");
hyA=read(file,"hy_0.r") + im*read(file,"hy_0.i");
close(file)

file=h5open("pj_fluxB.h5");
exB=read(file,"ex_0.r") + im*read(file,"ex_0.i");
eyB=read(file,"ey_0.r") + im*read(file,"ey_0.i");
hxB=read(file,"hx_0.r") + im*read(file,"hx_0.i");
hyB=read(file,"hy_0.r") + im*read(file,"hy_0.i");
close(file)

file=h5open("pj_fluxC.h5");
eyC=read(file,"ey_0.r") + im*read(file,"ey_0.i");
ezC=read(file,"ez_0.r") + im*read(file,"ez_0.i");
hyC=read(file,"hy_0.r") + im*read(file,"hy_0.i");
hzC=read(file,"hz_0.r") + im*read(file,"hz_0.i");
close(file)

file=h5open("pj_mode1.h5");
ex1=read(file,"ex.r") + im*read(file,"ex.i");
ey1=read(file,"ey.r") + im*read(file,"ey.i");
hx1=read(file,"hx.r") + im*read(file,"hx.i");
hy1=read(file,"hy.r") + im*read(file,"hy.i");
close(file)

file=h5open("pj_mode2.h5");
ex2=read(file,"ex.r") + im*read(file,"ex.i");
ey2=read(file,"ey.r") + im*read(file,"ey.i");
hx2=read(file,"hx.r") + im*read(file,"hx.i");
hy2=read(file,"hy.r") + im*read(file,"hy.i");
close(file)

file=h5open("pj_mode3.h5");
ex3=read(file,"ex.r") + im*read(file,"ex.i");
ey3=read(file,"ey.r") + im*read(file,"ey.i");
hx3=read(file,"hx.r") + im*read(file,"hx.i");
hy3=read(file,"hy.r") + im*read(file,"hy.i");
close(file)

file=h5open("pj_mode4.h5");
ex4=read(file,"ex.r") + im*read(file,"ex.i");
ey4=read(file,"ey.r") + im*read(file,"ey.i");
hx4=read(file,"hx.r") + im*read(file,"hx.i");
hy4=read(file,"hy.r") + im*read(file,"hy.i");
close(file)

sB=zeros(size(exB));
s1=zeros(size(ex1));
s2=zeros(size(ex2));
s3=zeros(size(ex3));
s4=zeros(size(ex4));
for n=1:length(sB)
  sB[n] = abs( conj(exB[n])*hyB[n] - conj(eyB[n])*hxB[n] );
  s1[n] = abs( conj(ex1[n])*hy1[n] - conj(ey1[n])*hx1[n] );
  s2[n] = abs( conj(ex2[n])*hy2[n] - conj(ey2[n])*hx2[n] );
  s3[n] = abs( conj(ex3[n])*hy3[n] - conj(ey3[n])*hx3[n] );
  s4[n] = abs( conj(ex4[n])*hy4[n] - conj(ey4[n])*hx4[n] );
end

clf();
subplot(6,4,1); plot(abs(exA[1,:,2])); 
subplot(6,4,2); plot(abs(eyA[1,:,2])); 
subplot(6,4,3); plot(abs(hxA[1,:,2])); 
subplot(6,4,4); plot(abs(hyA[1,:,2])); 

subplot(6,4,5); plot(abs(exB[1,:,2])); 
subplot(6,4,6); plot(abs(eyB[1,:,2])); 
subplot(6,4,7); plot(abs(hxB[1,:,2])); 
subplot(6,4,8); plot(abs(hyB[1,:,2])); 

subplot(6,4,9) ; plot(abs(ex1[1,:,2])); 
subplot(6,4,10); plot(abs(ey1[1,:,2])); 
subplot(6,4,11); plot(abs(hx1[1,:,2])); 
subplot(6,4,12); plot(abs(hy1[1,:,2])); 

subplot(6,4,13); plot(abs(ex2[1,:,2])); 
subplot(6,4,14); plot(abs(ey2[1,:,2])); 
subplot(6,4,15); plot(abs(hx2[1,:,2])); 
subplot(6,4,16); plot(abs(hy2[1,:,2])); 

subplot(6,4,17); plot(abs(ex3[1,:,2])); 
subplot(6,4,18); plot(abs(ey3[1,:,2])); 
subplot(6,4,19); plot(abs(hx3[1,:,2])); 
subplot(6,4,20); plot(abs(hy3[1,:,2])); 

subplot(6,4,21); plot(abs(ex4[1,:,2])); 
subplot(6,4,22); plot(abs(ey4[1,:,2])); 
subplot(6,4,23); plot(abs(hx4[1,:,2])); 
subplot(6,4,24); plot(abs(hy4[1,:,2])); 
