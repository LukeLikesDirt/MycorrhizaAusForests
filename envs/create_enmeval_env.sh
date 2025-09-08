#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

echo "Checking for existing enmeval environment..."
if conda info --envs | grep -q "enmeval"; then
    echo "Removing existing enmeval environment..."
    conda env remove -n enmeval -y
    echo "✓ Existing environment removed"
fi

echo "Creating conda environment from envs/enmeval.yml using mamba..."
mamba env create -f envs/enmeval.yml

if [ $? -eq 0 ]; then
    echo "✓ Environment created successfully!"
    source ~/.bashrc
    conda activate enmeval

    if [ $? -eq 0 ]; then
        echo "Installing R packages from GitHub/CRAN..."
        
        # ENMeval, DSM and DSMextra from CRAN/GitHub
        #Rscript -e "remotes::install_github('jamiemkass/ENMeval', ref = '2.0.4', dependencies=TRUE)"
        Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/ENMeval/ENMeval_2.0.4.tar.gz", repos = NULL, type = "source")'
        Rscript -e "install.packages('dsm', repos='https://cloud.r-project.org')"
        Rscript -e "remotes::install_github('densitymodelling/dsmextra', dependencies=TRUE)"

        echo "✓ All GitHub/CRAN packages installed successfully."
        echo "To use the environment in the future, run: conda activate enmeval"
    else
        echo "✗ Could not activate environment!"
    fi
else
    echo "✗ Error creating environment!"
    exit 1
fi
