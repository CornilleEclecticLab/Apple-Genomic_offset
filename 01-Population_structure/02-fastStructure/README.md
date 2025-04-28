# Run fastStructure

In jubail, create the environment with faststructure and the necessary dependencies
'''
module load miniconda
conda env create -f environment.yaml
'''

Create and submit faststructure file
'''
bash ./in/run_faststructure
'''
a bed file is asked to run faststructure, but used only the name without the extension