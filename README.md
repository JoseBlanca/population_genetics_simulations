# Population genetics simulations and practices

In this repository you'll find some practices related to population genetics that I have used in my teaching.

Each practice is located within a [Jupyter Notebook](https://jupyter.org/) file.

Practices implemented:

 - [Population definition using PCA](src/tomato_pop_definition_practice.ipynb)
 - [Genetic drift](src/genetic_drift_practice.ipynb)
 - [Three drifting populations](src/drifting_pops_practice.ipynb)
 - [Bottleneck](src/bottleneck_practice.ipynb)
 
 ## How to run the code

 To run the code you need [Jupyter](https://jupyter.org/) and some population genetic libraries that have been used to create the simulations:

  - [msprime](https://tskit.dev/msprime/docs/stable/intro.html)
  - [scikit-allel](https://scikit-allel.readthedocs.io/)
  - [demesdraw](https://github.com/grahamgower/demesdraw)

The full list of dependencies can be found in the files requirements.txt or environment.yml

If you are new to Python the easiest way to run the code in your computer is:

 - Install [Anaconda](https://www.anaconda.com/products/individual)
 - Launch "Anaconda Navigator"
 - Run a console inside Anaconda Navigator (in windows CMD.exe promt will do)
 - Inside the terminal run the command: 

```
conda install --yes -c conda-forge scikit-allel ipywidgets seaborn msprime plotly demesdraw jupyter-dash
```

 - Once the command is run close the terminal
 - Download and unzip the [practices](https://github.com/JoseBlanca/population_genetics_simulations/archive/refs/heads/main.zip)
 - In Anaconda Navigator open Launch "Jupyter Lab"
 - In the Jyputer Lab file browser go to the practices src folder and one any of the Jupyter Notebook files (the ones with the extension .ipynb)
 