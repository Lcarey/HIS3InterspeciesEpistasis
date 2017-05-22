import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg 
mpl.rcParams['pdf.fonttype'] = 42
import os
import numpy as np
from IPython.html.widgets.widget_float import FloatProgress
from IPython.display import display
import pandas as pd
from scipy import stats

def density_plot(x, y, nbins=42, log=False):
    mask = (~np.isnan(x)) & (~np.isnan(y))
    x = x[mask]
    y = y[mask]
    H, xedges, yedges = np.histogram2d(x,y,bins=nbins)
    ix = np.searchsorted(xedges, x)
    ix[ix == nbins] = nbins - 1
    iy = np.searchsorted(yedges, y)
    iy[iy == nbins] = nbins - 1
    v = H[ix, iy]
    i = v.argsort()
    cc = v[i]
    if log:
        cc = np.log(cc + 1)
    plt.scatter(x[i], y[i], c=cc, s=3, edgecolor='')


def plot_better(width=10, height=5, grid='xy', legend=False, visible_axes=True):
    plt.figure(figsize=(width, height))
    ax = plt.subplot(111)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",  
            labelbottom="on", left="off", right="off", labelleft="on")
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False) 
    if visible_axes:
        ax.spines["bottom"].set_visible(True) 
        ax.spines["bottom"].set_color('gray') 
        ax.spines["left"].set_visible(True)   
        ax.spines["left"].set_color('gray')
    else:
        ax.spines["bottom"].set_visible(False)  
        ax.spines["left"].set_visible(False) 
    
    if grid == 'xy':
        ax.xaxis.grid(True) 
        ax.yaxis.grid(True) 
    if grid == 'x':
        ax.xaxis.grid(True) 
    if grid == 'y':
        ax.yaxis.grid(True) 
    if legend:
        plt.legend(loc=2, bbox_to_anchor=(1.05, 1), frameon=False)

    return ax
    


def improve_plot(ax, grid='xy', legend=False, visible_axes=True):
    plt.tick_params(axis="both", which="both", bottom="off", top="off",  
            labelbottom="on", left="off", right="off", labelleft="on")
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False) 
    if visible_axes:
        ax.spines["bottom"].set_visible(True) 
        ax.spines["bottom"].set_color('gray') 
        ax.spines["left"].set_visible(True)   
        ax.spines["left"].set_color('gray')
    else:
        ax.spines["bottom"].set_visible(False)  
        ax.spines["left"].set_visible(False) 
    
    if grid == 'xy':
        ax.xaxis.grid(True) 
        ax.yaxis.grid(True) 
    if grid == 'x':
        ax.xaxis.grid(True) 
    if grid == 'y':
        ax.yaxis.grid(True) 
    if legend:
        plt.legend(loc=2, bbox_to_anchor=(1.05, 1), frameon=False)

    return ax


class Counter:

    def __init__(self, initial_count=0):
        self.count = initial_count


    def get_number(self, as_string=True):
        self.count += 1
        if as_string:
            return '%02d' % self.count
        return self.count

    def reset(self, initial_count=0):
        self.count = initial_count


def save_image(image_counter, title, folder, prefix):
    figure_name = '%s_img%s_%s.png' %(prefix, image_counter.get_number(), '_'.join(title.split()))
    plt.savefig(os.path.join(folder, figure_name), dpi=300)

def save_image_simply(title, folder):
    figure_name = '%s.png' %('_'.join(title.split()))
    plt.savefig(os.path.join(folder, figure_name), dpi=300)



def get_colors(number_of_colors, colormap):
    step = 255./(number_of_colors-1)
    return [colormap(int(i*step)) for i in range(number_of_colors)]


def hex_to_rgb(value):
    """Return (red, green, blue) for the color given as #rrggbb."""
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)




# # # # # # # # # # # # 

def get_AB(mut_combination, mutA, mutB):
    mut_combination = mut_combination.split(':')
    mut_combination.extend([mutA, mutB])
    mut_combination = list(set(mut_combination))
    mut_combination = sorted(mut_combination, key=lambda m: int(m[:-1]))
    return ':'.join(mut_combination)

def get_wt_mutA_mutB(mut_combination, mutA, mutB):
    mut_combination = mut_combination.split(':')
    wild_type_combination = ':'.join([m for m in mut_combination if m != mutA and m != mutB])
    mutA_combination = ':'.join([m for m in mut_combination if m != mutA])
    mutB_combination = ':'.join([m for m in mut_combination if m != mutB])
    return wild_type_combination, mutA_combination, mutB_combination


def get_foursomes(df_with_aaMutations_as_index, mutA, mutB):
    
    contains_A = find_genotypes_containing_mutations(df_with_aaMutations_as_index, [mutA])
    contains_B = find_genotypes_containing_mutations(df_with_aaMutations_as_index, [mutB])
    contains_AB = find_genotypes_containing_mutations(df_with_aaMutations_as_index, [mutA, mutB])

    if min([len(contains_A), len(contains_B), len(contains_AB)]) < 10:
        return pd.Panel()
    
    contains_A = contains_A[~contains_A.mut_list_Scer.apply(lambda s: contains_mutation(s, mutB))]
    contains_B = contains_B[~contains_B.mut_list_Scer.apply(lambda s: contains_mutation(s, mutA))]
    
    contains_A['wt'] = contains_A.mut_list_Scer.apply(lambda mut_combination: get_wt_mutA_mutB(mut_combination, mutA, mutB)[0])
    contains_B['wt'] = contains_B.mut_list_Scer.apply(lambda mut_combination: get_wt_mutA_mutB(mut_combination, mutA, mutB)[0])
    contains_AB['wt'] = contains_AB.mut_list_Scer.apply(lambda mut_combination: get_wt_mutA_mutB(mut_combination, mutA, mutB)[0])
    
    contains_A.set_index('wt', inplace=True)
    contains_B.set_index('wt', inplace=True)
    contains_AB.set_index('wt', inplace=True)

    foursome = pd.Panel.from_dict({'wild_type':df_with_aaMutations_as_index.set_index('wt'), 
                                   'mutA':contains_A, 'mutB':contains_B, 'mutAB':contains_AB}, intersect=True)
    return foursome


def get_foursomes_for_every_pair(df_with_aaMutations_as_index, mut_combinations, prefix, folder_to_save):
    f = FloatProgress(min=0, max=len(mut_combinations)+1)
    display(f)
    for mutA, mutB in mut_combinations:  
        foursome = get_foursomes(df_with_aaMutations_as_index, mutA, mutB)
        if len(foursome.major_axis) > 0:
            fn = prefix + 'foursome_mutA_%s_mutB_%s.hdf' %(mutA, mutB)
            foursome.to_hdf(os.path.join(folder_to_save, fn), 'data')
        f.value += 1


def remove_mutation_from_combination(mut_combination, mutation):
    return ':'.join([m for m in mut_combination.split(':') if m != mutation])



def check_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory



def get_mutation_in_all_backgrounds(df, mutation, lowest_acceptable_fitness=None):
    containing_mutation = find_genotypes_containing_mutations(df, mutation).copy()
    containing_mutation['wt'] = containing_mutation['mut_list_Scer'].apply(lambda mut_combination: remove_mutation_from_combination(mut_combination, mutation))
    containing_mutation.set_index('wt', inplace=True)

    if lowest_acceptable_fitness:
        df = df[df['s'] >= lowest_acceptable_fitness]
    panel = pd.Panel.from_dict({'wild_type':df.set_index('wt'), 'mutA':containing_mutation}, intersect=True)
    return panel

def get_fitness_impacts_in_all_backgrounds(df, mutation, lowest_acceptable_fitness=None):
    panel = get_mutation_in_all_backgrounds(df, mutation, lowest_acceptable_fitness=lowest_acceptable_fitness)
    return panel['mutA']['s'] - panel['wild_type']['s']


def consists_of_known_mutations(mutations, list_of_known_singles):
    for mutation in mutations.split(':'):
        if not mutation in list_of_known_singles:
            return False
    return True

def find_genotypes_containing_mutations(df, mutations):
    if type(mutations) == str:
        mutations = mutations.split(':')
    for m in mutations:
        df = df[df.mut_list_Scer.apply(lambda comb: contains_mutation(comb, m))]
    return df

def find_genotype(df, mutations_as_string):
    masked = df[df.mut_list_Scer == mutations_as_string]
    assert len(masked) < 2
    if len(masked) == 1:
        return masked
    else:
        return None
    
def contains_mutation(mutations, mutation):
    return mutation in mutations.split(':')




# # # # # # PyMol # # # # # # # #

import pymol
from pymol import cmd, stored
import matplotlib
from matplotlib import cm


def open_or_fetch(PDB_ID_or_filename, object_name=None):
    if len(PDB_ID_or_filename) in [4,5] and '.' not in PDB_ID_or_filename:
        cmd.fetch(PDB_ID_or_filename, async=0)
    else:
        if not object_name:
            object_name = PDB_ID_or_filename
        cmd.load(PDB_ID_or_filename, object_name)
        
def save_session(filename_pse='test.pse', pymol_viewer_version='1.72'):
    cmd.set('pse_export_version', pymol_viewer_version)
    cmd.save(filename_pse)

def white_and_beautiful(representation='cartoon'):
    cmd.hide('lines', 'all')
    cmd.show(representation, 'all')
    cmd.select('waters', 'name o')
    cmd.hide('everything', 'waters')
    cmd.color('white', 'all')
    cmd.set('bg_rgb', '(1,1,1)')
    cmd.set('surface_quality', '1')
    cmd.set('transparency', '0.5')
    cmd.set('ray_opaque_background', 'off')

    
# def prepare_GFP_2WUR():
#     cmd.fetch('2WUR', async=0)
#     white_and_beautiful()
#     cmd.select('waters', 'name o')
#     cmd.select('chr', 'resn GYS')
#     cmd.select('aa_64_68', 'resi 64+68')
#     cmd.select('aa_64_68_mainchain', 'aa_64_68 and name C+CO+CA+N')
#     cmd.hide('everything', 'waters')
#     cmd.show('sticks', 'chr')
#     cmd.color('green', 'chr')
#     cmd.show('sticks', 'aa_64_68_mainchain')
    
def color_positions(positions, values=None, representation='spheres', colormap=matplotlib.cm.cool, 
    constant_color=120, print_colors=False):
    colors_used = []
    # only positive values
    # if str(values) != 'None':
    #     assert min(values) >= 0
    if type(constant_color) == int or type(constant_color) == float:
        color = colormap(constant_color)
    elif type(constant_color) == str:
        color = mpl.colors.hex2color(constant_color)
    elif type(constant_color) == tuple:
        color = constant_color
    else:
        print 'Weird color!'

    if str(values) != 'None':
        values = np.array(values) - min(values)
        values = 1. * values / max(values)
    for index, position in enumerate(positions):
        if str(values) != 'None':
            color=colormap(values[index])
            colors_used.append((values[index], color))
        colorName = "color_" + str(position)
        selName = "temp_selection"
        cmd.set_color(colorName, color[0:3])
        cmd.select(selName, 'resi %s' %position)
        cmd.show(representation, selName)
        cmd.color(colorName, selName)
    if print_colors:
        return set(colors_used)
    else:
        return None
        
def get_residues_from_selection(selection_name, only_numbers=True):
    stored.list=[]
    if only_numbers:
        cmd.iterate("(%s & n. ca)" %selection_name, "stored.list.append(resi)")
        return [int(resi) for resi in stored.list]
    else:
        cmd.iterate("(%s & n. ca)" %selection_name, "stored.list.append((resi, resn))")
        return [(int(resi), resn) for resi, resn in stored.list]
    
def save_session_properly(session_counter, title, folder, prefix):
    session_name = '%s_pse%s_%s.pse' %(prefix, session_counter.get_number(), '_'.join(title.split()))
    save_session(os.path.join(folder, session_name))



# # # # # # # # # # # # 

def plot_segment_positions(ax, positions, segment_colors, scale=1):
    old_y = 2*scale
    for row in positions.iterrows():
        for position in row[1].positions_Uniprot_P06633:
            new_y = np.random.choice([1*scale, 2*scale])
            while new_y == old_y:
                new_y = np.random.choice([1*scale, 2*scale])
        x = row[1].positions_Uniprot_P06633
        plt.plot(x, [new_y for e in x], marker='o', markeredgecolor='none', lw=3, alpha=0.7, ls='none',
            label=row[1].segment, color=segment_colors[row[1].segment])
        # text_position = np.mean(max([v for v in x if v < np.median(x)]), min([v for v in x if v > np.median(x)]))
        text_position = np.median(x)
        plt.text(text_position, new_y+1, row[1].segment, va='center', ha='center')
        old_y = new_y


def printlist(l):
    for item in l:
        print item
    print


from pymol import cmd, stored

def interfaceResidues(cmpx, cA='c. A', cB='c. B', cutoff=1.0, selName="interface"):
    """
    interfaceResidues -- finds 'interface' residues between two chains in a complex.
    
    PARAMS
        cmpx
            The complex containing cA and cB
        
        cA
            The first chain in which we search for residues at an interface
            with cB
        
        cB
            The second chain in which we search for residues at an interface
            with cA
        
        cutoff
            The difference in area OVER which residues are considered
            interface residues.  Residues whose dASA from the complex to
            a single chain is greater than this cutoff are kept.  Zero
            keeps all residues.
            
        selName
            The name of the selection to return.
            
    RETURNS
        * A selection of interface residues is created and named
            depending on what you passed into selName
        * An array of values is returned where each value is:
            ( modelName, residueNumber, dASA )
            
    NOTES
        If you have two chains that are not from the same PDB that you want
        to complex together, use the create command like:
            create myComplex, pdb1WithChainA or pdb2withChainX
        then pass myComplex to this script like:
            interfaceResidues myComlpex, c. A, c. X
            
        This script calculates the area of the complex as a whole.  Then,
        it separates the two chains that you pass in through the arguments
        cA and cB, alone.  Once it has this, it calculates the difference
        and any residues ABOVE the cutoff are called interface residues.
            
    AUTHOR:
        Jason Vertrees, 2009.        
    """
    # Save user's settings, before setting dot_solvent
    oldDS = cmd.get("dot_solvent")
    cmd.set("dot_solvent", 1)
    
    # set some string names for temporary objects/selections
    tempC, selName1 = "tempComplex", selName+"1"
    chA, chB = "chA", "chB"
    
    # operate on a new object & turn off the original
    cmd.create(tempC, cmpx)
    cmd.disable(cmpx)
    
    # remove cruft and inrrelevant chains
    cmd.remove(tempC + " and not (polymer and (%s or %s))" % (cA, cB))
    
    # get the area of the complete complex
    cmd.get_area(tempC, load_b=1)
    # copy the areas from the loaded b to the q, field.
    cmd.alter(tempC, 'q=b')
    
    # extract the two chains and calc. the new area
    # note: the q fields are copied to the new objects
    # chA and chB
    cmd.extract(chA, tempC + " and (" + cA + ")")
    cmd.extract(chB, tempC + " and (" + cB + ")")
    cmd.get_area(chA, load_b=1)
    cmd.get_area(chB, load_b=1)
    
    # update the chain-only objects w/the difference
    cmd.alter( "%s or %s" % (chA,chB), "b=b-q" )
    
    # The calculations are done.  Now, all we need to
    # do is to determine which residues are over the cutoff
    # and save them.
    stored.r, rVal, seen = [], [], []
    cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')

    cmd.enable(cmpx)
    cmd.select(selName1, 'none')
    for (model,resi,diff) in stored.r:
        key=resi+"-"+model
        if abs(diff)>=float(cutoff):
            if key in seen: continue
            else: seen.append(key)
            rVal.append( (model,resi,diff) )
            # expand the selection here; I chose to iterate over stored.r instead of
            # creating one large selection b/c if there are too many residues PyMOL
            # might crash on a very large selection.  This is pretty much guaranteed
            # not to kill PyMOL; but, it might take a little longer to run.
            cmd.select( selName1, selName1 + " or (%s and i. %s)" % (model,resi))

    # this is how you transfer a selection to another object.
    cmd.select(selName, cmpx + " in " + selName1)
    # clean up after ourselves
    cmd.delete(selName1)
    cmd.delete(chA)
    cmd.delete(chB)
    cmd.delete(tempC)
    # show the selection
    cmd.enable(selName)
    
    # reset users settings
    cmd.set("dot_solvent", oldDS)
    
    return rVal

cmd.extend("interfaceResidues", interfaceResidues)

def flatten_list(nested_list):
    print 'Check this function'
    return [item for sublist in nested_list for item in sublist]


aa_1 = list("ACDEFGHIKLMNPQRSTVWYX")
aa_3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR XXX".split()
aa_long = ['alanine', 'cysteine', 'aspartic acid', 'glutamic acid', 'phenylalanine', 'glycine',
                'histidine', 'isoleucine', 'lysine', 'leucine', 'methionine', 'asparagine', 'proline',
                'glutamine', 'arginine', 'serine', 'threonine', 'valine', 'tryptophan', 'tyrosine',
                'unknown_aa']
aa_123 = dict(zip(aa_1, aa_3))
aa_321 = dict(zip(aa_3, aa_1))
