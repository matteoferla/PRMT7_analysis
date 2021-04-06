## M81del

> This is an exported Jupyter notebook

The model was made in a cluster-hosted Jupyter notebook with Pyrosetta and Remodel

### Load module

```python
import pyrosetta
from pyrosetta_help.init_ops import make_option_string, configure_logger

logger = configure_logger()
pyrosetta.distributed.maybe_init(extra_options=make_option_string(no_optH=False,
                                                ex1=None,
                                                ex2=None,
                                                #mute='all',
                                                ignore_unrecognized_res=True,
                                                load_PDB_components=False,
                                                ignore_waters=False)
                               )
```

### Load pose
```python
from pyrosetta_help.common_ops import pose_from_file

pose = pose_from_file('model2.pdb', params_filenames=['SAH.params'])

metal_setup = pyrosetta.rosetta.protocols.simple_moves.SetupMetalsMover()
metal_setup.set_remove_hydrogens(True)  # default
metal_sele = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(pyrosetta.rosetta.core.chemical.ResidueProperty.METAL)
metal_setup.set_metal_selector(metal_sele)
metal_setup.apply(pose)

scorefxn = pyrosetta.get_fa_scorefxn()
stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
stm.score_type_from_name('metalbinding_constraint')
scorefxn.set_weight(stm.score_type_from_name('metalbinding_constraint'), 1)
```

### Define blueprint
If Remodel fails to place the residues, it results in a model with the default amino acid (glycine in this case)
placed in structures matching if possible the requested secondary structure.
The `correct_and_relax` command below replaces the glycine in the quasi-C&alpha; trace and relaxes the neighbourhood.

```python
from pyrosetta_help.blueprint_maker import Blueprinter
blue = Blueprinter.from_pose(pose)
#blue.mutate(81, 'A')
#blue.mutate(r, 'A')
r = pose.pdb_info().pdb2pose(res=81, chain='A')
loopend = pose.pdb_info().pdb2pose(res=89, chain='A')

del blue[r]
blue[r-3:loopend] = 'NATAA'
blue.write('mut.blu')
blue.set('mut.blu')

blue.generic_aa = 'G'
blue.find_neighbors = True
remodel = blue.get_remodelmover(dr_cycles=5,
                               max_linear_chainbreak=0.2)


variant = pose.clone()
remodel.apply(variant)
blue.correct_and_relax(variant)
variant.dump_pdb('temp.pdb')
blue.show_poses_aligned(pose, variant)
```

<div style="font-family:monospace; display: inline-block; white-space: nowrap;">HQEIARSSYADMLHDKDRNVKYYQGIRAAVSRVKDRGQKALVLDIGTGTGLLSMMAVTAGADFCYAIEVFKPMADAAVKIVEKNGFSDKIKVINKHSTEVTVGPEGDMPCRANILVTELFDTELIGEGALPSYEHAHRHLVEENCEAVPHRATVYAQLVESGRMWSWNKLFPIHVQTSLGEQVIVPPVDVESCPGAPSVCDIQLNQVSPADFTVLSDVLPMFSIDFSKQVSSSAACHSRRFEPLTSGRAQVVLSWWDIEMDPEGKIKCTMAPFWAHSDPEEMQWRDHWMQCVYFLPQEEPVVQGSALYLVAHHDDYCVWYSLQRTSPEKNERVRQMRPVCDCQAHLLWNRPRFGEINDQDRTDRYVQALRTVLKPDSVCLCVSDGSLLSVLAHHLGVEQVFTVESSAASHKLLRKIFKANHLEDKINIIEKRPELLTNEDLQGRKVSLLLGEPFFTTSLLPWHNLYFWYVRTAVDQHLGPGAMVMPQAASLHAVVVEFRDLWRIRSPCGDCEGFDVHIMDDMIKRALDFRESREAEPHPLWEYPCRSLSEPWQILTFDFQQPVPLQPLCAEGTVELRRPGQSHAAVLWMEYHLTPECTLSTGLLEPADPEGGCCWNPHCKQAVYFFSPAPDPRALLGGPRTVSYAVEFHPDTGDIIMEFRHAZZ<br>.....................................................*..................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................<br>HQEIARSSYADMLHDKDRNVKYYQGIRAAVSRVKDRGQKALVLDIGTGTGLLS-MAVTAGADFCYAIEVFKPMADAAVKIVEKNGFSDKIKVINKHSTEVTVGPEGDMPCRANILVTELFDTELIGEGALPSYEHAHRHLVEENCEAVPHRATVYAQLVESGRMWSWNKLFPIHVQTSLGEQVIVPPVDVESCPGAPSVCDIQLNQVSPADFTVLSDVLPMFSIDFSKQVSSSAACHSRRFEPLTSGRAQVVLSWWDIEMDPEGKIKCTMAPFWAHSDPEEMQWRDHWMQCVYFLPQEEPVVQGSALYLVAHHDDYCVWYSLQRTSPEKNERVRQMRPVCDCQAHLLWNRPRFGEINDQDRTDRYVQALRTVLKPDSVCLCVSDGSLLSVLAHHLGVEQVFTVESSAASHKLLRKIFKANHLEDKINIIEKRPELLTNEDLQGRKVSLLLGEPFFTTSLLPWHNLYFWYVRTAVDQHLGPGAMVMPQAASLHAVVVEFRDLWRIRSPCGDCEGFDVHIMDDMIKRALDFRESREAEPHPLWEYPCRSLSEPWQILTFDFQQPVPLQPLCAEGTVELRRPGQSHAAVLWMEYHLTPECTLSTGLLEPADPEGGCCWNPHCKQAVYFFSPAPDPRALLGGPRTVSYAVEFHPDTGDIIMEFRHAZZ<br>  Score=663</div>

### Refine globally
```python
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
relax.apply(variant)
```

### Dump, score and check
```python
from pyrosetta_help.common_ops import add_bfactor_from_score
add_bfactor_from_score(variant)
scorefxn(variant) - scorefxn(pose) 
```


    27.281696571466
    

```python
import nglview as nv
view = nv.show_rosetta(variant)
view.clear_representations()
view.add_hyperball(selection='hetero')
view.add_tube(radiusType="bfactor", color="bfactor", radiusScale=0.01, colorScale="RdYlBu")
view
```

    NGLWidget()
    
None of the residues scores poorly individually, but the collective score of the neighbourhood goes up.

## Disregarded methods

The above approach assumes that the helix is maintained, but the span beyond 81 if rotated by 60˚.
However, in light of the bad score, the possibility that a hairpin is formed instead is investigated.


```python
# 'p.Met81del hairpin
def remove_termini(pose, preceding:int):
    # preceding precedes the discontinuity
    # LOWER_CONNECT N
    # UPPER_CONNECT C
    rm_upper = pyrosetta.rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue
    rm_lower = pyrosetta.rosetta.core.conformation.remove_lower_terminus_type_from_conformation_residue
    rm_upper(pose.conformation(), preceding)
    rm_lower(pose.conformation(), preceding)
    rm_lower(pose.conformation(), preceding + 1)
    rm_upper(pose.conformation(), preceding + 1)

def fix_loop(pose, start, stop, cutpoint):
    remove_termini(pose, cutpoint)
    # cutpoint: resi before discontinuity
    lm = pyrosetta.rosetta.protocols.loop_modeler.LoopModeler()
    loops = pyrosetta.rosetta.protocols.loops.Loops()
    loop = pyrosetta.rosetta.protocols.loops.Loop(start,
                                                  stop, 
                                                  cutpoint) #cutpoint
    #loop.auto_choose_cutpoint(pose)
    loops.add_loop(loop)
    lm.set_loops(loops)
    # these are enabled by default. here for quick changing.
    lm.enable_centroid_stage()
    lm.enable_fullatom_stage()
    lm.enable_build_stage()
    lm.apply(pose)
```

```python
variant = pose.clone()
r = variant.pdb_info().pdb2pose(res=81, chain='A')
variant.delete_residue_slow(r)
fix_loop(variant, r -3, r+2, r - 1)
```