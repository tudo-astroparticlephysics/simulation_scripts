import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
from icecube import dataclasses

### This function allows the user to plot a muon track ###

def muon_track_plotter(muon, tree, show=False, savefig=None):   
    # muon = I3Particle
    # tree = I3MCTree of muon
    # show = Boolean
    # savefig = String (name of plot)
    
    # Get direction and position data of muon
    dir_x, dir_y = muon.dir.x, muon.dir.y 
    pos_x, pos_y = muon.pos.x, muon.pos.y 

    # Get track and energies and types of daughter particles
    xs = []
    ys = []
    energies = []
    types = []

    for p in tree.get_daughters(muon):
        if p.type_string == 'PairProd' or p.type_string == 'DeltaE' or p.type_string == 'NuclInt' or p.type_string == 'Brems':
            xs.append(p.pos.x)
            ys.append(p.pos.y)
            energies.append(p.energy)
            types.append(p.type_string)

    # Define color map for different types of energy deposition
    color_map = {'Brems':"#FF6600",
                 'PairProd':"#00FF33",
                 'NuclInt':"#FFFF33",
                 'DeltaE':"#330066"}

    # Plot energy deposition 
    fig, ax = plt.subplots()

    max_r = 20
    max_val = max(np.log10(energies))

    for i in range(len(energies)):
        c = plt.Circle((xs[i], ys[i]), radius=np.log10(energies[i])/max_val*max_r, color=color_map[types[i]])
        fig.gca().add_artist(c)
    
    ax.set_aspect('equal')
    plt.xlabel('x / m')
    plt.ylabel('y / m')

    # Plot muon track
    track_xs = np.linspace(xs[0], xs[-1], 100)
    track_ys = (track_xs - pos_x) / dir_x * dir_y + pos_y 
    plt.plot(track_xs, track_ys, '-', color='red')

    # Plot hadronic cascade from the prom. interaction  
    cascade = tree.get_daughters(muon)[1]
    c = plt.Circle((cascade.pos.x,cascade.pos.y),radius=np.log10(cascade.energy)/max_val*max_r,color="magenta")
    fig.gca().add_artist(c)

    # Create legend with mpatches
    brems_patch = mpatches.Patch(color=color_map['Brems'], label='Brems')
    pair_patch = mpatches.Patch(color=color_map['PairProd'], label='PairProd')
    nuclint_patch = mpatches.Patch(color=color_map['NuclInt'], label='NuclInt')
    deltae_patch = mpatches.Patch(color=color_map['DeltaE'], label='DeltaE')
    cascade_patch = mpatches.Patch(color="magenta", label='Cascade')

    plt.legend([brems_patch,pair_patch,nuclint_patch,deltae_patch, cascade_patch],["Brems","PairProd","NuclInt","DeltaE", "Cascade"],
           prop={"size":8},bbox_to_anchor=(1, 0.5),loc="center left")



    # Save figure
    if savefig != None:
        os.system('mkdir -p plots')
        plt.savefig('plots/{}.pdf'.format(savefig))

    # Show figure
    if show == True:
        plt.show()

### This functions allows to plot the track with a tray module ###

def muon_plotter(frame):
    muon = dataclasses.get_most_energetic_muon(frame['I3MCTree'])
    muon_track_plotter(muon, frame['I3MCTree'], show=False, savefig='muon_track{}'.format(muon.minor_id))


### This function searches for the biggest energy deposition ###

def max_energy_depo_new_tree(muon, tree, new_tree):
    # muon = I3Particle
    # tree = I3MCTree of muon

    # Get daughters of muon
    daughters = tree.get_daughters(muon)

    # Get energy depositions
    energies = []

    for d in daughters:
        energies.append(d.energy)

    # Get index of energy list with maximum energy
    max_index = energies.index(max(energies))

    # Get particle in daughters with maximum energy
    max_particle = daughters[max_index]

    print('-----This is the particle with the maximum energy:----- \n {}'.format(max_particle))
    print('-----This is the old tree:----- \n {}'.format(tree))


    # Add primary to new tree
    new_tree.add_primary(max_particle)


    print('-----This is the new tree:----- \n {}'.format(new_tree))

    # Get access to max particle via tree 
    print(tree.get_particle(max_particle))


### Function max energy depo new tree for tray ###

def max_energy_new_tree(frame):
    frame.Put('NewTree', dataclasses.I3MCTree())
    muon = dataclasses.get_most_energetic_muon(frame['I3MCTree'])
    max_energy_depo_new_tree(muon, frame['I3MCTree'], frame['NewTree'])

### Define functions for muon splitting ###


def muon_splitter(tray, name):
    from icecube import dataclasses, icetray

    def dir_pos_tree(frame):
        frame.Put('NewDirection', dataclasses.I3Direction())
        frame.Put('NewPosition', dataclasses.I3Position())		
        frame.Put('NewTree', dataclasses.I3MCTree())


    tray.AddModule(dir_pos_tree, 'addDirPos', Streams=[icetray.I3Frame.DAQ])
