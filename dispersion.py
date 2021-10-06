#!/usr/bin/env python

###############################################################
#                                                             #
#                D I S P E R S I O N . P Y                    #
#                                                             #
###############################################################
''' 
       ALTERNATIVE CODE FOR PLOTTING BANDSTRUCTURES 
                 FROM A CASTEP .BANDS FILE
'''

# Let us import all the stuff we need, shouldnt require any specialist packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from fractions import Fraction
import sys
import os
from itertools import cycle
import argparse
import ase.io as io
import ase.dft.bz as bz
import warnings

def blockPrint():
    sys.stdout = open(os.devnull, 'w')
  
# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

# Define some constants
hartree = 27.211386245988
fracs=np.array([0.5,0.0,0.25,0.75,0.33333333,0.66666667])


# pdos reader
def pdos_read(seed,species):
    from scipy.io import FortranFile as FF

    f=FF(seed+'.pdos_bin', 'r','>u4')
    
    version=f.read_reals('>f8')
    header=f.read_record('a80')[0]
    num_kpoints=f.read_ints('>u4')[0]
    num_spins=f.read_ints('>u4')[0]
    num_popn_orb=f.read_ints('>u4')[0]
    max_eigenvalues=f.read_ints('>u4')[0]
    
    orbital_species=f.read_ints('>u4')
    orbital_ion=f.read_ints('>u4')
    orbital_l=f.read_ints('>u4')
    
    kpoints=np.zeros((num_kpoints,3))
    pdos_weights=np.zeros((num_popn_orb,max_eigenvalues,num_kpoints,num_spins))
    for nk in range(0,num_kpoints):
        record=f.read_record('>i4','>3f8')
        kpt_index,kpoints[nk,:]=record
        for ns in range(0,num_spins):
            spin_index=f.read_ints('>u4')[0]
            num_eigenvalues=f.read_ints('>u4')[0]

            for nb in range(0,num_eigenvalues):
                pdos_weights[0:num_popn_orb,nb,nk,ns]=f.read_reals('>f8')
                
                #norm=np.sqrt(np.sum((pdos_weights[0:num_popn_orb,nb,nk,ns])**2))
                norm=np.sum((pdos_weights[0:num_popn_orb,nb,nk,ns]))
                pdos_weights[0:num_popn_orb,nb,nk,ns]=pdos_weights[0:num_popn_orb,nb,nk,ns]/norm

    if species:
        num_species=len(np.unique(orbital_species))
        pdos_weights_sum=np.zeros((num_species,max_eigenvalues,num_kpoints,num_spins))

        for i in range(0,num_species):
            loc=np.where(orbital_species==i+1)[0]
            pdos_weights_sum[i,:,:,:]=np.sum(pdos_weights[loc,:,:,:],axis=0)
                    
            
    else:
        num_orbitals=4
        pdos_weights_sum=np.zeros((num_orbitals,max_eigenvalues,num_kpoints,num_spins))
        pdos_colours=np.zeros((3,max_eigenvalues,num_kpoints,num_spins))

        r=np.array([1,0,0])
        g=np.array([0,1,0])
        b=np.array([0,0,1])
        k=np.array([0,0,0])
                                
        
        
        for i in range(0,num_orbitals):
            loc=np.where(orbital_l==i)[0]
            if len(loc)>0:

                pdos_weights_sum[i,:,:,:]=np.sum(pdos_weights[loc,:,:,:],axis=0)
        #print(kpoints[1])
        #for nb in range(num_eigenvalues):
        #    print(pdos_weights_sum[:,nb,1,0])
    pdos_weights_sum=np.where(pdos_weights_sum>1,1,pdos_weights_sum)
    pdos_weights_sum=np.where(pdos_weights_sum<0,0,pdos_weights_sum)
    return np.round(pdos_weights_sum,7)


def path_finder():
    
    # Open the cell
    path_str=bv_latt.special_path

    path_points=[]
    path_labels=[]
    for L in path_str:
        if L==",":
            break
        path_labels.append(L)
        path_points.append(special_points[L])

    print("%BLOCK SPECTRAL_KPOINT_PATH")
    for i in range(len(path_labels)):
        print("%.5f %.5f %.5f" %(path_points[i][0],path_points[i][1],path_points[i][2]),"#",path_labels[i])
    print("%ENDBLOCK SPECTRAL_KPOINT_PATH")
    
    

    

def cart_to_abc(lattice):
     a=np.sqrt( lattice[0,0]**2+lattice[0,1]**2+lattice[0,2]**2)
     b=np.sqrt( lattice[1,0]**2+lattice[1,1]**2+lattice[1,2]**2)
     
     c=np.sqrt( lattice[2,0]**2+lattice[2,1]**2+lattice[2,2]**2)
     alpha=( lattice[1,0]* lattice[2,0]+lattice[1,1]* lattice[2,1]+lattice[1,2]* lattice[2,2])/(b*c)
     alpha=np.arccos(alpha)
     beta =( lattice[2,0]* lattice[0,0]+lattice[2,1]* lattice[0,1]+lattice[2,2]* lattice[0,2])/(c*a)
     beta =np.arccos(beta)
     gamma=( lattice[0,0]* lattice[1,0]+lattice[0,1]* lattice[1,1]+ lattice[0,2]*lattice[1,2])/(a*b)
     gamma=np.arccos(gamma)
     return a,b,c,alpha,beta,gamma






def calc_phonons(buff_seed):
    no_ions = 0
    no_kpoints = 0
    no_branches = 0
    no_electrons = 0
    unit = 0 

    # Open the phonon file
    phonon_file=buff_seed+".phonon"
    phonon=open(phonon_file,'r')

    lines=phonon.readlines()

    no_ions=int(lines[1].split()[-1])
    no_branches=int(lines[2].split()[-1])
    no_kpoints=int(lines[3].split()[-1])


    lattice=np.zeros((3,3))
    lattice[0]=[i for i in lines[8].split()]
    lattice[1]=[i for i in lines[9].split()]
    lattice[2]=[i for i in lines[10].split()]

    
    #make the arrays
    energy_array=np.empty(shape=(no_kpoints,no_branches))

    kpoint_array=np.empty(shape=(no_kpoints)) # the array holding the number of the kpoint
    kpoint_list=[] # array of the kpoint vectors


    kpoint_string=lines[15::no_branches+3+no_ions*no_branches]
    
        
    for i in range(len(kpoint_string)):
        kpoint_array[i]=int(kpoint_string[i].split()[1])
        
        #Empty list for vectors
        vec=[]
        vec.append(float(kpoint_string[i].split()[2]))
        vec.append(float(kpoint_string[i].split()[3]))
        vec.append(float(kpoint_string[i].split()[4]))
        kpoint_list.append(vec)
        
        #   print(vec)
        #Lets get the eigen values into the big array

    for k in range(0,no_kpoints):
        
        ind=16 + (k) * (3+no_branches+no_ions*no_branches)


        
        
        
        energy_array[k,:]=np.array([float(i.split()[-1]) for i in lines[ind:ind+no_branches]])
    

    sort_array=kpoint_array.argsort()
    
    kpoint_list=np.array(kpoint_list)[sort_array]


    return energy_array,sort_array,kpoint_list,kpoint_array,no_kpoints,no_ions,lattice


    

 
    
# Variables we need from the bands file

def calc_bands(buff_seed):
    no_spins = 0
    no_kpoints = 0
    fermi_energy = 0
    no_electrons = 0
    no_electrons_2 = 0
    no_eigen = 0
    no_eigen_2 = 0

    # Open the bands file
    bands_file=buff_seed+".bands"
    bands=open(bands_file,'r')

    lines=bands.readlines()

    no_spins=int(lines[1].split()[-1])
    no_kpoints=int(lines[0].split()[-1])
    fermi_energy=float(lines[4].split()[-1])

    if no_spins==1:
        fermi_energy=float(lines[4].split()[-1])
        no_electrons =float(lines[2].split()[-1])
        no_eigen  = int(lines[3].split()[-1])

    if no_spins==2:
        spin_polarised=True
        no_electrons=float(lines[2].split()[-2])
        no_electrons_2=float(lines[2].split()[-1])
        no_eigen  = int(lines[3].split()[-2])
        no_eigen_2=int(lines[3].split()[-1])

    lattice=np.zeros((3,3))
    lattice[0]=[i for i in lines[6].split()]
    lattice[1]=[i for i in lines[7].split()]
    lattice[2]=[i for i in lines[8].split()]
    lattice=lattice/1.889
    
    #make the arrays
    energy_array=np.empty(shape=(no_kpoints,no_eigen))
    energy_array_2=np.empty(shape=(no_kpoints,no_eigen_2))

    kpoint_array=np.empty(shape=(no_kpoints)) # the array holding the number of the kpoint
    kpoint_list=[] # array of the kpoint vectors

    if no_spins==1:
        kpoint_string=lines[9::no_eigen+2]
    else:
        kpoint_string=lines[9::no_eigen+3+no_eigen_2]
    
        #loop through the kpoints to split it
        
    for i in range(len(kpoint_string)):
        kpoint_array[i]=int(kpoint_string[i].split()[1])
        
        #Empty list for vectors
        vec=[]
        vec.append(float(kpoint_string[i].split()[2]))
        vec.append(float(kpoint_string[i].split()[3]))
        vec.append(float(kpoint_string[i].split()[4]))
        kpoint_list.append(vec)
        
        #   print(vec)
        #Lets get the eigen values into the big array

    for k in range(0,no_kpoints):
        if no_spins==1:
            ind=9+k*no_eigen+2*(k+1)
            energy_array[k,:]=hartree*np.array([float(i)-fermi_energy for i in lines[ind:ind+no_eigen]])
        if no_spins==2:
            ind=9+k*(no_eigen+no_eigen_2+1)+2*(k+1)
        

            energy_array[k,:]=hartree*np.array([float(i)-fermi_energy for i in lines[ind:ind+no_eigen]])
            energy_array_2[k,:]=hartree*np.array([float(i)-fermi_energy for i in lines[ind+no_eigen+1:ind+no_eigen+1+no_eigen_2]])

    sort_array=kpoint_array.argsort()

    kpoint_list=np.array(kpoint_list)[sort_array]


    return energy_array,energy_array_2,sort_array,kpoint_list,kpoint_array,no_spins,no_kpoints,fermi_energy,no_electrons,no_electrons_2,no_eigen,no_eigen_2,lattice




def check_sym(vec):
    frac=[]
    for i in vec:
        #frac.append(i.as_integer_ratio()[0])
        #frac.append(i.as_integer_ratio()[1])
        buff=[]
        for j in fracs:
            buff.append(np.isclose(i,j))
        frac.append(any(buff))
    
        
    
    if all(frac):
        #print(vec)
        return True
    else:
        return False 


def main():
    warnings.filterwarnings("ignore")
    #matplotlib.rcParams['mathtext.fontset'] = 'stix'
    #matplotlib.rcParams['font.family'] = 'STIXGeneral'
    #matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
    #matplotlib.use('macOsX')
    matplotlib.rc('text', usetex = True)
    plt.style.use("classic")
    
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    
    
    #Do the parser
    parser = argparse.ArgumentParser(description= "Utillity for plotting bandstructurs from a CASTEP run.")
    parser.add_argument("seed",help="The seed from the CASTEP calculation.")
    parser.add_argument("--save",action="store_true",help="Save DOS as .pdf with name <seed>-dos.pdf.")
    parser.add_argument("-m","--multi",action="store_true",help="Set lines multicoloured.")
    parser.add_argument("-l","--line",help="Set linewidth.",default=0.75)
    parser.add_argument("--lim",help="Provide plotting limits around the Fermi energy.",nargs=2,default=[None,None])
    parser.add_argument("-s","--spin",help="Plot spin-up and spin-down channels.",action="store_true")
    parser.add_argument("-d","--debug",action='store_true',help="Debug flag.")
    #parser.add_argument("--sym",help="Provide crystal symmetry for plot labels.",default=None)
    parser.add_argument("--soc",help="Seedname of second bands file containing SOC bandstructure.",default=None)
    parser.add_argument("--n_up",help="Indices of up bands to be highlighted",nargs="+")
    parser.add_argument("--n_down",help="Indices of down bands to be highlighted",nargs="+")
    parser.add_argument("-f","--flip",action="store_true",help="Plot with a global spin flip")
    parser.add_argument("--fontsize",help="Font size",default=20)
    parser.add_argument("--title",help="Add a title for saving")
    parser.add_argument("--fig",help="add figure caption")
    parser.add_argument("-e","--exe",help="File extension for saving",default="png")
    parser.add_argument("--dos",help="Prodide some data files for DOS plots adjoining bandstructure",nargs="+")
    parser.add_argument("--path",help="Compute a suitable band path for the cell and exit.",nargs="*")
    parser.add_argument("--pdos",help="Use .pdos_bin file to project orbital information",action='store_true')
    parser.add_argument("--species",help="Project pdos onto species rather than orbitals",action='store_true')
    parser.add_argument("--phonon",help="Plot phonon dispersion curve",action='store_true')
    parser.add_argument("-b","--bandgap",help="Indicate bandgap on plots",action="store_true")
    parser.add_argument("--no_plot",help="Supress plotting of dispersions",action="store_true")
    args = parser.parse_args()
    seed = args.seed
    save = args.save
    multi= args.multi
    linewidth=np.float(args.line)
    lim= args.lim
    debug=args.debug
    spin_split=args.spin
    #sym=args.sym
    SOC=args.soc
    spin_polarised=False
    n_up=args.n_up
    n_down=args.n_down
    flip=args.flip
    text=float(args.fontsize)
    title=args.title
    fig_cap=args.fig
    exe=args.exe
    dos_files=args.dos
    path=args.path
    pdos=args.pdos
    species=args.species
    do_phonons=args.phonon
    bg=args.bandgap
    no_plot=args.no_plot
    blockPrint()

    print(do_phonons)
    
    # Dothe path and labels
    cell=io.read(seed+".cell")
    bv_latt=cell.cell.get_bravais_lattice()
    special_points=bv_latt.get_special_points()
    atoms=np.unique(cell.get_chemical_symbols())[::-1]
    
    enablePrint()


    
    if path==[]:
        path_finder()
        sys.exit()
    else:
        if path!=None:
            path_points=[]
            path_labels=[]
            
            for i in path:
                
                try:
                    path_point=special_points[i]
                except:
                    print()
                    print("Error: %s has no symmetry point %s"%(bv_latt.name,i))
                    sys.exit()
                    path_points.append(path_point)
                    path_labels.append(i)

            print("%BLOCK SPECTRAL_KPOINT_PATH")
            for j in range(len(path_labels)):
                print("%.5f %.5f %.5f" %(path_points[j][0],path_points[j][1],path_points[j][2]),"#",path_labels[j])
            print("%ENDBLOCK SPECTRAL_KPOINT_PATH")
            sys.exit()
        
            
    if n_up!=None:
        n_up=np.array(n_up,dtype=int)-1
        spin_split=False
    else:
        n_up=[]
    if n_down!=None:
        n_down=np.array(n_down,dtype=int)-1
        spin_split=False
    else:
        n_down=[]
    if SOC != None:
        doSOC=True
    else :
        doSOC=False
    if dos_files!=None:
        do_dos=True
    else:
        do_dos=False


    
    bands_file=True
    
    if multi and spin_split:
        multi=False
        
    if doSOC:
        multi=False
        spin_split=False
    
    #set the colours
    if spin_split:
        spin_up="r"
        spin_do="b"
    if flip:
        spin_up="b"
        spin_do="r"
    
    else :
        spin_up="black"
        spin_do="black"


    #calculate the pdos if needed
    if pdos:
        pdos_weights=pdos_read(seed,species)
        
    if doSOC:
        energy_array_soc,energy_array_2,sort_array_soc,kpoint_list_soc,kpoint_array_soc,no_spins,no_kpoints,fermi_energy,no_electrons,no_electrons_2,no_eigen,no_eigen_2,lattice2=calc_bands(SOC)

    if not do_phonons:
        energy_array,energy_array_2,sort_array,kpoint_list,kpoint_array,no_spins,no_kpoints,fermi_energy,no_electrons,no_electrons_2,no_eigen,no_eigen_2,lattice=calc_bands(seed)
    
        if energy_array_2.shape[1]!=0:


            vb_max_up=np.max(energy_array[:,int(no_electrons)-1])
            vb_max_down=np.max(energy_array_2[:,int(no_electrons_2)-1])
            cb_min_up=np.min(energy_array[:,int(no_electrons)])
            cb_min_down=np.min(energy_array_2[:,int(no_electrons_2)])
            
            band_gap_up=cb_min_up-vb_max_up
            band_gap_down=cb_min_down-vb_max_down
            print("Band gap (up)   : %6.3f eV"%band_gap_up)
            print("Band gap (down) : %6.3f eV"%band_gap_down)
            vb_max_ind_up=np.where(energy_array[sort_array][:,int(no_electrons)-1]==vb_max_up)[0][-1]
            vb_max_ind_down=np.where(energy_array_2[sort_array][:,int(no_electrons_2)-1]==vb_max_down)[0][-1]
            
            cb_min_ind_up=np.where(energy_array[sort_array][:,int(no_electrons)]==cb_min_up)[0][-1]
            cb_min_ind_down=np.where(energy_array_2[sort_array][:,int(no_electrons_2)]==cb_min_down)[0][-1]
            
            k_max_loc_up=kpoint_array[sort_array][vb_max_ind_up]
            k_max_loc_down=kpoint_array[sort_array][vb_max_ind_down]
            k_min_loc_up=kpoint_array[sort_array][cb_min_ind_up]
            k_min_loc_down=kpoint_array[sort_array][cb_min_ind_down]
            
            
        else:
            vb_max=np.max(energy_array[:,int(no_electrons/2)-1])
            cb_min=np.min(energy_array[:,int(no_electrons/2)])
            
            band_gap=cb_min-vb_max
            print("Band gap : %6.3f eV"%band_gap)
            vb_max_ind=np.where(energy_array[sort_array][:,int(no_electrons/2)-1]==vb_max)[0][-1]
            cb_min_ind=np.where(energy_array[sort_array][:,int(no_electrons/2)]==cb_min)[0][-1]
            
            k_max_loc=kpoint_array[sort_array][vb_max_ind]
            k_min_loc=kpoint_array[sort_array][cb_min_ind]
            
    else:
        energy_array,sort_array,kpoint_list,kpoint_array,no_kpoints,no_ions,lattice=calc_phonons(seed)

    
    
    a,b,c,alpha,beta,gamma=cart_to_abc(lattice)
    a1,a2,a3=lattice[0],lattice[1],lattice[2]
    b1=2*np.pi*np.cross(a2,a3)/(np.dot(a1,np.cross(a2,a3)))
    b2=2*np.pi*np.cross(a3,a1)/(np.dot(a1,np.cross(a2,a3)))
    b3=2*np.pi*np.cross(a1,a2)/(np.dot(a1,np.cross(a2,a3)))
    kalpha=np.arccos(np.dot(a2,a3)/(np.linalg.norm(a2)*np.linalg.norm(a3)))
    kbeta=np.arccos(np.dot(a1,a3)/(np.linalg.norm(a1)*np.linalg.norm(a3)))
    kgamma=np.arccos(np.dot(a2,a1)/(np.linalg.norm(a2)*np.linalg.norm(a1)))
    
    #matplotlib.rc('text', usetex = True)
    
    # Here we do the analysis of the kpoints and the symmetry.. It's going to be horific!
    
    #define all the greek letters we will use for weird ones

    
    
    if no_plot:
        sys.exit()
    k_ticks=[]
    for i,vec in enumerate(kpoint_list):
        if check_sym(vec):
            k_ticks.append(kpoint_array[i])
    
    tol=1e-5
    tol=[tol,tol,tol]
    
    kpoint_grad=[]
    for i in range(1,len(kpoint_list)):
        diff=kpoint_list[i]-kpoint_list[i-1]
        kpoint_grad.append(diff)
    
    kpoint_2grad=[]
    high_sym=[0]
    for i in range(1,len(kpoint_grad)):
        diff=kpoint_grad[i]-kpoint_grad[i-1]
        kpoint_2grad.append(diff)
        #print(diff)
    
        if any(np.abs(diff)>tol):
    
            # print(diff)
            high_sym.append(i)
    high_sym.append(len(kpoint_list)-1)
    high_sym=np.array(high_sym)+1
    
    ##################### SOC ###################
    if doSOC:
        k_ticks_soc=[]
        for i,vec in enumerate(kpoint_list_soc):
            if check_sym(vec):
                k_ticks_soc.append(kpoint_array_soc[i])
    
        tol=1e-5
        tol=[tol,tol,tol]
    
        kpoint_grad_soc=[]
        for i in range(1,len(kpoint_list_soc)):
            diff=kpoint_list_soc[i]-kpoint_list_soc[i-1]
            kpoint_grad_soc.append(diff)
    
        kpoint_2grad_soc=[]
        high_sym_soc=[0]
        for i in range(1,len(kpoint_grad_soc)):
            diff=kpoint_grad_soc[i]-kpoint_grad_soc[i-1]
            kpoint_2grad_soc.append(diff)
            #print(diff)
    
            if any(np.abs(diff)>tol):
    
            # print(diff)
                high_sym_soc.append(i)
        high_sym_soc.append(len(kpoint_list_soc)-1)
        high_sym_soc=np.array(high_sym_soc)+1
    
        #############################################
    
        if len(high_sym)!=len(high_sym_soc):
            print("SOC Bandsstructure Does not match")
            sys.exit()
    
        for i in range(1,len(high_sym)):
            high_up=int(high_sym[i])
            high_low=int(high_sym[i-1])
            soc_up=int(high_sym_soc[i])
            soc_low=int(high_sym_soc[i-1])
    
            nsoc=len(kpoint_array_soc[soc_low:soc_up])+1
            nhigh=len(kpoint_array[high_low:high_up])+1
            kpoint_array_soc[soc_low-1:soc_up]=np.linspace(high_low,high_up,nsoc,endpoint=True)
    
    
                    
    
    
    
    # Set up the plotting environment
    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif',weight='bold')
    
    #Do the fonts
    #matplotlib.rcParams['font.sans-serif'] = "Times New Roman"#Comic Sans MS"
    # Then, "ALWAYS use sans-serif fonts"
    #matplotlib.rcParams['font.family'] = "sans-serif"
    
    if not do_dos:
        fig, ax = plt.subplots(figsize=(7,7))
    else:
        from matplotlib.ticker import MaxNLocator
        fig, (ax, ax2) = plt.subplots(1, 2,sharey=True, gridspec_kw={'hspace': 0,'wspace': 0,'width_ratios': [2.4, 1]},figsize=(11,7))
        for file in dos_files:
            pdos_dat=np.loadtxt(file)
            shape=pdos_dat.shape[1]
            
            energy = pdos_dat[:,0]
            if lim[0]!= None:
                mask = (energy >= float(lim[0])) & (energy <= float(lim[1]))
            else:
                ax2.set_ylim(lim[0],lim[1])
                mask=[True]*len(energy)
                [mask]
    
            if shape==3:
                ax2.plot(pdos_dat[:,1][mask],energy[mask],linewidth=linewidth,linestyle="--",color="black")
            if shape==5:
                ax2.plot(2*(pdos_dat[:,1][mask]-pdos_dat[:,2][mask]),energy[mask],linewidth=linewidth,color="black")
    
        ax2.axhline(0,color="0.6",dashes=[8, 8],linewidth=1,)
        ax2.tick_params(axis='both', which='major', labelsize=text,length=7)
        ax2.set_xlabel(r"$\mathit{g}(\mathit{E}$) (states/eV)",fontsize=text)
        ax2.xaxis.set_major_locator(MaxNLocator(4)) 
        dos_ticks=ax2.get_xticks()
        dos_ticks=np.delete(dos_ticks,0)
        ax2.set_xticks(dos_ticks)
    
            
    for vline in high_sym:
        ax.axvline(vline,color="black",linewidth=1)
    ax.set_xticks(high_sym)
    ax.axhline(0,color="0.6",dashes=[8, 8],linewidth=1,)
    if not do_phonons:
        ax.set_ylabel(r'$\mathit{E}$-$\mathit{E}_{\mathrm{F}}$ (eV)',fontsize=text)
    else:
        ax.set_ylabel(r'$\omega$ (cm$^{-1}$)',fontsize=text)
    ax.set_xlim(1,no_kpoints)
    ax.tick_params(axis='both', which='major', labelsize=text,length=7)
    
    if lim[0]!= None:
        ax.set_ylim(float(lim[0]),float(lim[1]))
    
    
    
    #set the x labels
    ticks= []
    tol=1e-4
    '''
    if sym==None:
        for vec in kpoint_list[high_sym-1]:
            ticks.append("("+str(Fraction(vec[0]).limit_denominator())+","+str(Fraction(vec[1]).limit_denominator())+","+str(Fraction(vec[2]).limit_denominator())+")")
    
    
        ax.set_xticklabels(ticks)
        for tick in ax.get_xticklabels():
            tick.set_rotation(-30)'''
    
    ticks=[""]*len(high_sym)
    found=False
    for k_count,k in enumerate(kpoint_list[high_sym-1]):
        found=False
        for i in special_points:#sym_dict[sym]:
            
            #if abs(sym_dict[sym][i][0]-k[0])<tol and abs(sym_dict[sym][i][1]-k[1])<tol and abs(sym_dict[sym][i][2]-k[2])<tol:
            if abs(special_points[i][0]-k[0])<tol and abs(special_points[i][1]-k[1])<tol and abs(special_points[i][2]-k[2])<tol:
                if i=="G":
                    ticks[k_count]="$\Gamma$"
                else:
                    ticks[k_count]=i
                found=True
                #if not found:
                #    ticks.append("")
    
    ax.set_xticklabels(ticks)
    
    
        
    #plt.gcf().subplots_adjust(bottom=0.2)
    
    n_colors=cycle(['blue','red','green','black','purple','orange','yellow','cyan'])
    
    if bg:
        if energy_array_2.shape[1]!=0:
            ax.plot([k_max_loc_up,k_max_loc_up],[vb_max_up,cb_min_up],color='r',linewidth=linewidth*2)
            ax.plot([k_max_loc_down,k_max_loc_down],[vb_max_down,cb_min_down],color='b',linewidth=linewidth*2)
                    
            ax.plot([k_max_loc_up,k_min_loc_up],[cb_min_up,cb_min_up],color='r',linewidth=linewidth*2)
            ax.plot([k_max_loc_down,k_min_loc_down],[cb_min_down,cb_min_down],color='b',linewidth=linewidth*2)
            ax.text(k_max_loc_up*1.05,vb_max_up+(-vb_max_up+cb_min_up)*0.8/2,"%4.2f eV"%band_gap_up,fontsize=text)
            ax.text(k_max_loc_down*1.05,vb_max_down+(-vb_max_down+cb_min_down)*0.8/2,"%4.2f eV"%band_gap_down,fontsize=text)
    
    
        else:
            #ax.scatter(k_min_loc,cb_min)
            #ax.scatter(k_max_loc,vb_max)
            ax.plot([k_max_loc,k_max_loc],[vb_max,cb_min],color='k',linewidth=linewidth*2)
            ax.plot([k_max_loc,k_min_loc],[cb_min,cb_min],color='k',linewidth=linewidth*2)
            ax.text(k_max_loc*1.05,vb_max+(-vb_max+cb_min)*0.8/2,"%4.2f eV"%band_gap,fontsize=text)
    
    if multi:
        if not do_phonons:
            ax.plot(kpoint_array[sort_array],energy_array[sort_array],linewidth=linewidth)
            if no_spins==2:
                
                ax.plot(kpoint_array[sort_array],energy_array_2[sort_array])
        else:
            ax.plot(kpoint_array[sort_array],energy_array[sort_array],linewidth=linewidth)
    elif not do_phonons:
        if pdos:
            from matplotlib import colors
            from matplotlib.colors import ListedColormap
            from matplotlib.lines import Line2D
            import matplotlib.collections as mcoll
            import matplotlib.path as mpath
            
            def colorline(
                    x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0),
                    linewidth=3, alpha=1.0):
                """
                http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
                http://matplotlib.org/examples/pylab_examples/multicolored_line.html
                Plot a colored line with coordinates x and y
                Optionally specify colors in the array z
                Optionally specify a colormap, a norm function and a line width
                """
                
                # Default colors equally spaced on [0,1]:
                if z is None:
                    z = np.linspace(0.0, 1.0, len(x))
                    
                z = np.asarray(z)
                    
                segments = make_segments(x, y)
                
                lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                                          linewidth=linewidth, alpha=alpha)
                
                ax.add_collection(lc)
                
                return lc
            
            
            def make_segments(x, y):
                """
                Create list of line segments from x and y coordinates, in the correct format
                for LineCollection: an array of the form numlines x (points per line) x 2 (x
                and y) array
                """
                
                points = np.array([x, y]).T.reshape(-1, 1, 2)
                
                
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                
                return segments
            
            
    
            if species:
                n_cat=len(atoms)
            else:
                n_cat=4
                
                
            basis=[]
            for i in range(n_cat):
                basis.append(np.array(colors.to_rgba(next(n_colors))))
    
        
        
            for nb in range(no_eigen):               
                # calculate the colour
                cmap_array=np.zeros((len(kpoint_array),4))
                for i in range(n_cat):
                    
                    cmap_array[:,0]+=pdos_weights[i,nb,:,0]*basis[i][0]#/n_cat
                    cmap_array[:,1]+=pdos_weights[i,nb,:,0]*basis[i][1]#/n_cat
                    cmap_array[:,2]+=pdos_weights[i,nb,:,0]*basis[i][2]#/n_cat
                    cmap_array[:,3]+=pdos_weights[i,nb,:,0]*basis[i][3]#/n_cat
                    
                    #cmap_array[:,0:3]=cmap_array[:,0:3]/n_cat
                    cmap_array=np.where(cmap_array>1,1,cmap_array)
                    cmap = ListedColormap(cmap_array)
                
                z = np.linspace(0, 1, len(kpoint_array))
            
                colorline(kpoint_array[sort_array], energy_array[sort_array][:,nb], z, cmap=cmap, linewidth=3)
                ax.plot(kpoint_array[sort_array],energy_array[sort_array][:,nb],linewidth=linewidth,alpha=0)
            
                if no_spins==2:
        
                    for nb in range(no_eigen):               
                        # calculate the colour
                        cmap_array=np.zeros((len(kpoint_array),4))
                        for i in range(n_cat):
                            
                            cmap_array[:,0]+=pdos_weights[i,nb,:,1]*basis[i][0]#/n_cat
                            cmap_array[:,1]+=pdos_weights[i,nb,:,1]*basis[i][1]#/n_cat
                            cmap_array[:,2]+=pdos_weights[i,nb,:,1]*basis[i][2]#/n_cat
                            cmap_array[:,3]+=pdos_weights[i,nb,:,1]*basis[i][3]#/n_cat
                            
                            #cmap_array[:,0:3]=cmap_array[:,0:3]/n_cat
                        cmap_array=np.where(cmap_array>1,1,cmap_array)
                        cmap = ListedColormap(cmap_array)
                    
                    z = np.linspace(0, 1, len(kpoint_array))
            
                    colorline(kpoint_array[sort_array], energy_array_2[sort_array][:,nb], z, cmap=cmap, linewidth=3)
                    ax.plot(kpoint_array[sort_array],energy_array[sort_array][:,nb],linewidth=linewidth,alpha=0)
    
        
            custom_lines = []
            labels=[]
            for i in range(n_cat):
                custom_lines.append(Line2D([0], [0], color=basis[i], lw=3))
                if species:
                    labels.append(atoms[i])
                else:
                    labels=["s","p","d","f"]
                
                        
                    #custom_lines = [Line2D([0], [0], color=cmap(0.), lw=4),
                    #                Line2D([0], [0], color=cmap(.5), lw=4),
                    #                Line2D([0], [0], color=cmap(1.), lw=4)]
        
        
            ax.legend(custom_lines,labels,fontsize=text)
        
            
        else:
            ax.plot(kpoint_array[sort_array],energy_array[sort_array],color=spin_up,label="without SOC",linewidth=linewidth)
            for i in n_up:
                ax.plot(kpoint_array[sort_array],energy_array[sort_array][:,i],linewidth=linewidth,color=next(n_colors))
            c=1
            if no_spins==2:
                ax.plot(kpoint_array[sort_array],energy_array_2[sort_array],color=spin_do,label="without SOC",linewidth=linewidth)
                for i in n_down:
                    ax.plot(kpoint_array[sort_array],energy_array_2[sort_array][:,i],linewidth=linewidth,color=next(n_colors))
            
    
            if doSOC:
                #kpoint_array_soc=1+(kpoint_array[-1]-1)*(kpoint_array_soc-1)/(kpoint_array_soc[-1]-1)
                ax.plot(kpoint_array_soc,energy_array_soc[sort_array_soc],color=spin_do,label="with SOC",linewidth=linewidth,linestyle="--")
                handles, labels = plt.gca().get_legend_handles_labels()
                by_label = dict(zip(labels, handles))
                if not do_dos:
                    plt.legend(by_label.values(), by_label.keys(),loc="upper right",fontsize=text)
    
    
    
    else: #This is the part where we plot the phonons
        ax.plot(kpoint_array[sort_array],energy_array[sort_array],color=spin_up,label="without SOC",linewidth=linewidth)
      
                    
        
    
    if spin_polarised and debug:
        split_en=np.mean(energy_array-energy_array_2,axis=0)
            
    
    if not do_dos:
        plt.figtext(0.95, 0.96, fig_cap, wrap=True, horizontalalignment='center', fontsize=text)
    else:
    
        x=ax2.get_xlim()[1]*0.9
        y=ax2.get_ylim()[1]*0.85
        ax2.text(x,y,fig_cap,wrap=True, horizontalalignment='center', fontsize=text)
    title_seed=seed#.replace("_","\_")
    if save:
        if title!=None:
            plt.suptitle(title,fontsize=text)
    
        if do_phonons:
            plt.tight_layout()
            fig.savefig(seed+"-phonon."+exe)
        elif doSOC:
            plt.tight_layout()
            fig.savefig(seed+"-SOC-bs."+exe)
        elif do_dos:
            plt.tight_layout()
            fig.savefig(seed+"-SOC-bs-dos."+exe)
        else:
            plt.tight_layout()
            fig.savefig(seed+"-bs."+exe)
    else:
        plt.title(title_seed,fontsize=20)
        plt.tight_layout()
        plt.show()


if __name__=='__main__':
    main()
