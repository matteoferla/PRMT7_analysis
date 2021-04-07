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

```python
from pyrosetta_help.common_ops import pose_from_file

#pose = pose_from_file('model2.pdb', params_filenames=['SAH.params'])
pose = pose_from_file('model_cart_nov16.pdb', params_filenames=['SAH.params'])
metal_setup = pyrosetta.rosetta.protocols.simple_moves.SetupMetalsMover()
metal_setup.set_remove_hydrogens(True)  # default
metal_sele = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(pyrosetta.rosetta.core.chemical.ResidueProperty.METAL)
metal_setup.set_metal_selector(metal_sele)
metal_setup.apply(pose)

# ref2015: 
# scorefxn = pyrosetta.get_fa_scorefxn()
# parker2016:
pyrosetta.rosetta.basic.options.set_boolean_option('corrections:beta_nov16', True)
scorefxn = pyrosetta.create_score_function('beta_nov16_cart')

stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
scorefxn.set_weight(stm.score_type_from_name('metalbinding_constraint'), 1)
```

```python
from pyrosetta_help.score_mutants import MutantScorer
import os

MutantScorer.output_folder = 'nov16'
if not os.path.exists(MutantScorer.output_folder):
    os.mkdir(MutantScorer.output_folder)

model = MutantScorer(pose, modelname='model')
model.scorefxn = scorefxn 
model.strict_about_starting_residue = True

mutations = ['p.Trp494Arg', 'p.Arg32Thr',
'p.Arg387Gly', 'p.Arg497Gln',
'p.Glu94Lys' , 
 'p.E125G',
'p.Cys571Arg'] # 'p.Met81del' manually.

raw = model.score_mutations(mutations,
                            chain='A',
                            interfaces=(), #
                            preminimise=False,
                            distance=12,
                            cycles=5)
```

```python
import pandas as pd
from pyrosetta_help.score_mutants import extend_scores

scores = pd.DataFrame(raw)
extend_scores(scores)
scores.to_csv('scores_opt16_cart.csv')
```