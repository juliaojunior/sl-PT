import meep as mp
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math


# all lenghts in units of microns (um).
# PC periodicity is 1 um


def main(args):
    resolution = 20      # pixels/um
    eps = 13             # epsilon of waveguide
    cols = 3            # metade da quantidade de colunas no cristal
    lines  = 1  #7 é valor anterior            # metade da quantidade de linhas no cristal
    w = 1.2              # width of the  waveguide
    r = 0.2             # radius of holes
    N = args.N           # number of holes on either side of defect
    #pad = 1              # padding between last hole and PML
    dpml = 1             # PML thickness

    sy = 2*(dpml + lines)       # tamanho da celula em Y (perpend ao guia)
    #sy = args.sy        # size of cell in Y (perpend to wvg)
    fcen = args.fcen     # pulse centger frequency
    df = args.df         # pulse frequency width

    #sx = 2*(pad + dpml ) + lines      #  size of cell in X
    sx = 2*(dpml + cols)

    cell = mp.Vector3(sx,sy,0)


    blk = mp.Block(size=mp.Vector3(mp.inf,mp.inf,mp.inf), material=mp.Medium(epsilon=eps))
    geometry = [blk]

    for i in range(lines):
        geometry.append(mp.Cylinder(r,center=mp.Vector3(0,i)))
        for j in range(cols):
               geometry.append(mp.Cylinder(r,center=mp.Vector3(i,j)))   #      j,i)))
               geometry.append(mp.Cylinder(r,center=mp.Vector3(-i,j)))  #      j,-i)))
               geometry.append(mp.Cylinder(r,center=mp.Vector3(i,-j)))  #      -j,i)))
               geometry.append(mp.Cylinder(r,center=mp.Vector3(-i,-j))) #      -j,-i)))

    geometry.append(mp.Block(size=mp.Vector3(mp.inf,w,mp.inf),material=mp.Medium(epsilon=eps)))


    pml_layers = [mp.PML(1.0)]


    src = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                     component=mp.Ey,
                     center=mp.Vector3(-0.5*sx+dpml), # fonte na esquerda; para colocar na direita usar 0.5*sx - dpml
                     size=mp.Vector3(0,w))]


    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        sources=src,                 #symmetries=sym,
                        resolution=resolution)

    freg = mp.FluxRegion(center=mp.Vector3(0.5*sx-dpml-0.5),  # fluxo na direita; para colocar na esquerda usar -0.5*sx+dpml+0.5
                         size=mp.Vector3(0,2*w))
    nfreq = 500       # number of frequencies computed fluxes

    trans = sim.add_flux(fcen, df, nfreq, freg)     # transmitted flux



    vol = mp.Volume(mp.Vector3(0), size=mp.Vector3(sx))
    volTime = mp.Volume(mp.Vector3(((sx)/2 - dpml),0), size=mp.Vector3(0,1,0))   # dados para pulso no tempo no final do guia

    hvals = []

    def gethvals(sim):
        hvals.append(sim.get_array(center=mp.Vector3(), size=cell, component=mp.Hz))




    sim.run(mp.at_beginning(mp.output_epsilon),
            mp.in_volume(volTime,mp.to_appended("TimeVsE2",mp.at_every(0.4,mp.output_dpwr))),
            mp.in_volume(vol,mp.to_appended("hz-slice",mp.at_every(0.4,gethvals,mp.output_hfield_z))),
            until=300)


    sim.display_fluxes(trans)    # print flux spectrum  - linha realocada, ver abaixo

    epsilon0 = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
    plt.figure()
    plt.imshow(epsilon0.transpose(), interpolation='spline36', cmap='RdBu')
    plt.axis('off')
    plt.savefig('epsilon.png',format='png')
    plt.show()

#    sim.display_fluxes(trans)    # print flux spectrum     - posicao original da linha,  ver acima




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-N', type = int, default=3, help='number of holes on either side of defect')
    parser.add_argument('-Fator', type = int, default=1, help='multiplos da parte imaginaria do epsilon')
    parser.add_argument('-sy', type = int, default=6, help='size of cell in y direction (perpendicular to wvg)')
    parser.add_argument('-fcen', type = float, default=0.25, help='pulse center frequency')
    parser.add_argument('-df', type = float, default = 0.2, help = 'pulse frequency width')
    args = parser.parse_args()
    main(args)
