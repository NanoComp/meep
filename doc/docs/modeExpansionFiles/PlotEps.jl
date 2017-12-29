##################################################
# julia routine to produce a 3D plot of a .h5 file
# intended for use in visualizing 3D MEEP geometries
##################################################
using HDF5
using PyPlot

function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})
  m, n = length(vy), length(vx)
  vx = reshape(vx, 1, n)
  vy = reshape(vy, m, 1)
  (repmat(vx, m, 1), repmat(vy, 1, n))
end

function PlotEps(;filename= "eps-000000.00.h5", dy=1, alpha=0.2,
                 slice=false,
                 XMin=-4, XMax=4, YMin=-4, YMax=4, ZMin=-8, ZMax=8)

  file = h5open(filename);
  eps  = read(file, "eps");
  close(file);

  NX=size(eps,3);
  NY=size(eps,2);
  NZ=size(eps,1);

  clf()

  if (slice)
    imshow(eps[:,Int(NY/2),:],cmap="coolwarm");
    return
  end

  X=collect(linspace(XMin,XMax,NX));
  Y=collect(linspace(YMin,YMax,NY));
  Z=collect(linspace(ZMin,ZMax,NZ));
  (X,Z)=meshgrid(X,Z);

  subplot(111,projection="3d");
  for ny=1:dy:NY
   contourf(X,Z,eps[:,ny,:],zdir="z",offset=Y[ny],cmap="coolwarm", alpha=alpha)
   xlabel(L"$x$",fontsize=40)
   ylabel(L"$z$",fontsize=40)
   zlabel(L"$y$",fontsize=40)
   zlim([-4,4])
  end
end
