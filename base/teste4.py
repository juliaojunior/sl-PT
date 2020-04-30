#  Programa para testar funcao 
# onde o epsilon e' definido e tem
# valor imaginário. A funcao e' escrita e o dieletrico
# da simulacao e' baseada nela
# SEGUNDA VERSAO
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
    w = 1.2              # width of the  waveguide
    dpml = 1             # PML thickness
    largC = 16           # largura da celula
    altC = 16            # altura da celula

    sx = largC + dpml
    sy = altC + dpml

    fcen = args.fcen     # pulse centger frequency
    df = args.df         # pulse frequency width


    cell = mp.Vector3(sx,sy,0)


    def epsP(p):
       valorIm = -0.5 + math.pow(math.cos((p.y)*(2*math.pi/(altC/4))),2)
       return mp.Medium(epsilon=3, D_conductivity=2*math.pi*fcen*(valorIm)/3)
    epsP.do_averaging = True



    blk = mp.Block(size=mp.Vector3(mp.inf,mp.inf,mp.inf), material=mp.Medium(epsilon=eps))
    geometry = [blk]



#    AQUI ESTÁ A LINHA DE TESTE: USAMOS A FUNCAO PARA DEFINIR O EPSILON
#    geometry.append(mp.Block(size=mp.Vector3(mp.inf,w,mp.inf)))
    geometry.append(mp.Block(center=mp.Vector3(),size=mp.Vector3(mp.inf,altC,mp.inf),material=epsP))




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


    sim.run(mp.at_beginning(mp.output_epsilon),
            mp.to_appended("hz",mp.at_every(0.4,mp.output_hfield_z)),
            until=300)



    epsilon0 = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
    plt.figure()
    plt.imshow(epsilon0.transpose(), interpolation='spline36', cmap='RdBu')
    plt.axis('off')
    plt.savefig('epsilon.png',format='png')
    plt.show()





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fcen', type = float, default=0.25, help='pulse center frequency')
    parser.add_argument('-df', type = float, default = 0.2, help = 'pulse frequency width')
    args = parser.parse_args()
    main(args)
