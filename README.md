# BCrystal
BCrystal: An Interpretable Sequence-Based Protein Crystallization Predictor

## Motivation
Protein crystallization allows to study molecular structure. Novel in silico, accurate, sequence-based protein crystallization predictors are highly sought.

# Installation

### Requirements

This step will install all the dependencies required for running BCrystal. You do not need sudo permissions for this step.

  - Install Anaconda
    1. Download Anaconda (64 bit) installer python3.x for linux : https://www.anaconda.com/distribution/#download-section
    2. Run the installer : `bash Anaconda3-2019.03-Linux-x86_64.sh` and follow the instructions to install.
    3. Install xgboost: conda install -c conda-forge xgboost 
    4. Install shap: conda install -c conda-forge shap 
    5. Install Bio: conda install -c anaconda biopython 

  - R requirements
    - Run R REPL by running the following: `R`
    -  Install R libraries
       1.  Interpol (do `install.packages('Interpol')` )
       2.  bio3d    (do `install.packages('bio3d')` )
       3.  doParallel (do `install.packages('doParallel')`)
       4.  zoo      (do `install.packages('zoo')`)
       
    Quit R REPL: `quit()` 
 
  - SCRATCH (version SCRATCH-1D release 1.2) (http://scratch.proteomics.ics.uci.edu, Downloads: http://download.igb.uci.edu/#sspro)
    1. Run `wget http://download.igb.uci.edu/SCRATCH-1D_1.2.tar.gz`
    2. Run `tar -xvzf SCRATCH-1D_1.2.tar.gz`
    3. Run `cd SCRATCH-1D_1.2`
    4. Run `perl install.pl`
    5. Run `cd ..`
    6. Replace the blast in `SCRATCH-1D_1.2/pkg/blast-2.2.26` with a 64 bit version of `blast-2.2.26` if you are running on a 64 bit machine (`ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.26/`).
    
  - DISOPRED (version 3.16) (http://bioinfadmin.cs.ucl.ac.uk/downloads/DISOPRED/)
    1. Run `wget http://bioinfadmin.cs.ucl.ac.uk/downloads/DISOPRED/DISOPRED3.16.tar.gz`
    2. Run `tar -xvzf DISOPRED3.16.tar.gz`
    3. Run `cd DISOPRED/src/`
    4. Run `make clean; make; make install`

DISOPRED and SCRATCH-1D_1.2 should be in the same directory as Data folder. Data folder has the training and the 3 test set proteins in fasta format as well as files corresponding to their true labels - crystallized (1) or not (0).

# Run BCrsytal on New Test file

To run BCrystal on your own protein sequences you need the following three things:

  1. Protein Sequence File: Protein sequence/sequences of interest in fasta format (https://en.wikipedia.org/wiki/FASTA_format). We provide `data/Seq_solo.fasta` as an example 
  2. SCRATCH: Software used to extract structural features from a given protein sequence file. Follow instructions in the previous section to install SCRATCH.
  3. DISOPRED: Software used to extract disorder features from a given protein sequence file. Follow instructions in the previosu section to install DISOPRED.
  

### Execute in the command line
 
  1. `Rscript --vanilla features_PaRSnIP_v2.R <your-test>.fasta`
  2. `python xgb.py features.csv <your-test>.fasta <output_folder>`
  
  
## In the <output_folder> you will find 2 outputs
  1. prediction.csv - Containing the crystallization propensity
  2. bar_plot_i.png - where i=1 if a solo sequence is passed in fasta otherwise its the nth sequence in the test fasta file.
  
