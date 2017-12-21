using HDF5
using PyPlot
using PyCall
@pyimport matplotlib.animation as anim

##################################################
##################################################
##################################################
function plotImage(n, files, cmpt, extent)
  try
    file = h5open( String(files[n]) );
    data = read(file,string(cmpt,".r")) + im*read(file,string(cmpt,".i"))
    close(file); 
    imshow(abs(data[end:-1:1,:]),extent=extent)
    xlabel(L"$x$",fontsize=40)
    ylabel(L"$z$",fontsize=40)
    draw();
  catch
  end
end

##################################################
##################################################
##################################################
#function init(n,files, cmpt, extent)
#  return plotImage(1, files, cmpt, extent)
#end

function animate(n, files, cmpt, extent)
  return plotImage(n, files, cmpt, extent) 
end

##################################################
# fileList is a text file containing the name of one h5 file per line
##################################################
function makeMovie(h5fileList; cmpt="hz",
                               extent=[-6, 6, -8, 8],
                               interval=150,
                               mp4File="hz.mp4")

  pygui(true);
  files=readdlm(h5fileList);
  plotImage(1,files,cmpt,extent);

  fig = figure();
  movie = anim.FuncAnimation(fig, 
                             n->plotImage(n,files,cmpt,extent),
                             frames=length(files),
                             interval=interval, 
                             repeat=false);

  movie[:save](mp4File)
end
