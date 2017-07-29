# Python port of meep/examples/bend-flux.ctl
# From the Meep tutorial: transmission around a 90-degree waveguide bend in 2d.
from __future__ import division

import meep as mp
import meep.geom as gm
from meep.source import GaussianSource


# dummy material function needed to pass to structure()
# constructor as a placeholder before we can call
# set_materials_from_geometry
def dummy_eps(vec):
    return 1.0


def bend_flux(no_bend):

    sx = 16.0        # size of cell in X direction
    sy = 32.0        # size of cell in Y direction
    pad = 4.0        # padding distance between waveguide and cell edge
    w = 1.0          # width of waveguide
    resolution = 10  # (set-param! resolution 10)

    gv = mp.voltwo(sx, sy, resolution)
    gv.center_origin()
    the_structure = mp.structure(gv, dummy_eps, mp.pml(1.0))

    wvg_ycen = -0.5 * (sy - w - 2.0 * pad)  # y center of horiz. wvg
    wvg_xcen = 0.5 * (sx - w - 2.0 * pad)  # x center of vert. wvg

    e1 = gm.Vector3(1.0, 0.0, 0.0)
    e2 = gm.Vector3(0.0, 1.0, 0.0)
    e3 = gm.Vector3(0.0, 0.0, 1.0)

    dielectric = gm.Medium(epsilon_diag=gm.Vector3(12, 12, 12))
    if no_bend:
        center = gm.Vector3(y=wvg_ycen)
        size = gm.Vector3(float('inf'), w, float('inf'))
        objects = [gm.Block(size, e1, e2, e3, material=dielectric, center=center)]
        mp.set_materials_from_geometry(the_structure, objects)
    else:
        objects = []
        center = gm.Vector3(-0.5 * pad, wvg_ycen)
        size = gm.Vector3(sx - pad, w, float('inf'))
        objects.append(gm.Block(size, e1, e2, e3, material=dielectric, center=center))

        center = gm.Vector3(wvg_xcen, 0.5 * pad)
        size = gm.Vector3(w, sy - pad, float('inf'))
        objects.append(gm.Block(size, e1, e2, e3, material=dielectric, center=center))
        mp.set_materials_from_geometry(the_structure, objects)

    f = mp.fields(the_structure)

    fcen = 0.15  # pulse center frequency
    df = 0.1
    src = GaussianSource(fcen, df)
    v = mp.volume(mp.vec(1.0 - 0.5 * sx, wvg_ycen), mp.vec(0.0, w))
    f.add_volume_source(mp.Ez, src.swigobj, v)

    f_start = fcen - 0.5 * df
    f_end = fcen + 0.5 * df
    nfreq = 100  # number of frequencies at which to compute flux

    if no_bend:
        trans_volume = mp.volume(mp.vec(0.5 * sx - 1.5, wvg_ycen), mp.vec(0.0, 2.0 * w))
    else:
        trans_volume = mp.volume(mp.vec(wvg_xcen, 0.5 * sy - 1.5), mp.vec(2.0 * w, 0.0))

    trans_vl = mp.volume_list(trans_volume, mp.Sz)
    trans = f.add_dft_flux(trans_vl, f_start, f_end, nfreq)

    refl_volume = mp.volume(mp.vec(-0.5 * sx + 1.5, wvg_ycen), mp.vec(0.0, 2.0 * w))
    refl_vl = mp.volume_list(refl_volume, mp.Sz)
    refl = f.add_dft_flux(refl_vl, f_start, f_end, nfreq)

    dataname = "refl-flux"
    if not no_bend:
        refl.load_hdf5(f, dataname)
        refl.scale_dfts(-1.0)

    eval_point = mp.vec(0.5 * sx - 1.5, wvg_ycen) if no_bend else mp.vec(wvg_xcen, 0.5 * sy - 1.5)
    deltaT = 50.0
    next_check_time = f.round_time() + deltaT
    tol = 1.0e-3
    max_abs = 0.0
    cur_max = 0.0
    done = False

    while not done:
        f.step()

        # manually check fields-decayed condition
        absEz = abs(f.get_field(mp.Ez, eval_point))
        cur_max = max(cur_max, absEz)
        if f.round_time() >= next_check_time:
            next_check_time += deltaT
            max_abs = max(max_abs, cur_max)
            if max_abs > 0.0 and cur_max < tol * max_abs:
                done = True
            cur_max = 0.0

        # printf("%.2e %.2e %.2e %.2e\n",f.round_time(),absEz,max_abs,cur_max)

    if no_bend:
        refl.save_hdf5(f, dataname)

    print("{}\t\t | {}\t\t | {}".format("Time", "trans flux", "refl flux"))
    f0 = fcen - 0.5 * df
    fstep = df / (nfreq - 1)
    trans_flux = trans.flux()
    refl_flux = refl.flux()
    for nf in range(nfreq):
        print("{}\t\t | {}\t\t | {}".format(f0 + nf * fstep, trans_flux[nf], refl_flux[nf]))


def main(args):

    bend_flux(True)
    bend_flux(False)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
