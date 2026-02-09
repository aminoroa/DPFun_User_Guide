# Intro
Following is a guideline on how to use DPFunc (https://github.com/CSUBioGroup/DPFunc) on HPCC of Michigan State University. 
For reference to DPFunc and its applications, please read the realted paper: https://www.nature.com/articles/s41467-024-54816-8

We demonstrate the instructions to predict the function of proteins using the already-developed DPFunc model. To achieve this, we will 
    1. Clone DPFunc from github, 
    2. Provide the model and model data
    3. Prepare the input files for 7 different protein sequences 
    4. Run the model for the these 7 protein sequences
    5. visualize and save the results

# 1. Clone DPFunc from github

## 1.1. connect to HPCC through your terminal 

SSH <username@hpcc.msu.edu>
Enter password

## 1.2. Load a proper development node
** we recommend using dev-amd20-v100

SSH load dev-amd20-v100

## 1.3. Clone DPFunc from GitHub in the home directory

cd ~
git clone https://github.com/CSUBioGroup/DPFunc.git


## 1.4. Create the Conda environment necessary to run DPFunc
module purge
module load Miniforge3
module load CUDA/12.6.0
cd DPFunc
conda env create -f DPFunc_env.yml
conda activate py38dgl

# 2. Provide the model and model data

# 2.1. providing the model

Download the trained model from this link: https://drive.google.com/file/d/1V0VTFTiB29ilbAIOZn0okBQWPlbOI3wN/view?usp=drive_link

Upload this file to ~/DPFunc/ Directory
tar -xvf save_models.tar
** executing the command above should created ~/DPFunc/save_models directory with 9 .pt files related to bp, cc, and mf modules of the model

# 2.2. Provide the model data

Download the model data from this link: https://drive.google.com/file/d/1qrxbkk450GJzhVfqnAN9Ms798owrq96n/view?usp=sharing​

Upload the dataset in ./LargeScaleData directory

# 3. Prepare the input files for 7 different protein sequences

Here we demonstrate how you should prepare the required input files for your protein of interests. We would need to create input files below so that we can run the model.

3.1. PDB structure files
3.2. Residue-level ESM features
3.3. InterPro domain features
3.4. Graph features 

As a demonstration, we did this for 7 proteins mentioned below:
-3H3B
-3KDM
-4MN8
-8HND
-8IKW
-8IQS
-8JYR

## 3.1. PDB structure files

-From protein data bank (https://www.rcsb.org), download the the file related to each protein with "Legacy PDB" format. Place the data within ~/DPFunc/data/PDB/PDB_folder directory.
So,  ~/DPFunc/data/PDB/PDB_folder directory should contain 7 <protein_name>.pdb files.

-Then, we need to created .pkl files from .pdb files.
-To do so, add create_test_pid_list.py to ~/DPFunc directory.
-then, run the script:
python create_test_pid_list.py

-After successfully running the command above, you will have test_pid_list.pkl in ~/DPFunc/data/ directory.

# 3.2. Residue-level ESM features

-use make_pdb_seqs_for_myproteins.py

-run the script above which genenates pdb_seqs.pkl in ~/DPFunc/processed_file directory.

python make_pdb_seqs_for_myproteins.py

-then use process_esm.py
** you will likly need to install fair-esm (pip install fair-esm)

-run the script above which generates esm data: est_part_0.pkl in ~/DPFunc/processed_file/esm_emds directory:

python process_esm.py
 

## 3.3. InterPro Domain features 

DPFunc employ interpro to extranct the domains related to each protein.
If you do:

module avail 2>&1 | grep -i interpro

you will find out that Interpro scan with various version are already available on HPCC:
    InterProScan/5.62-94.0-foss-2022b 
    InterProScan/5.72-103.0-foss-2023b (D) 
    InterProScan_data/5.62-94.0-foss-2022b 
    InterProScan_data/5.72-103.0-foss-2023b

-first add convert_to_fasta.py to ~/DPFunc directory

-then execute this script to convert pdb_seq.pkl files to .fasta files.

-then submit run_interproscan_myproteins.sb to HPCC in order to get the Interproscan TSV output.

-Then, you will need to convert interproscan tsv file (~/DPFunc/processed_data/interproscan_out/myproteins.interproscan.tsv) into DPFunc interpro vectors.

-to do this use iprscan_tsv_to_dpfunc_vectors.py

-run the iprscan_tsv_to_dpfunc_vectors.py:

python iprscan_tsv_to_dpfunc_vectors.py

** running the script above successfully will generate 7 <protein name>.pkl file at ~/DPFunc/processed_file/interpro_22369 directory

-run make_pdb_points_for_myproteins.py 

python make_pdb_points_for_myproteins.py 

-run make_ca_aligned_seqs_and_points to align pdb_seq.pkl and pdb_point.pkl files

python make_ca_aligned_seqs_and_points

-rerun process_esm.py for the corrected sequence

python process_esm.py

** check the misalignment using misalignment_check.py. If you get "misalignment=[]" it means everything is fine to proceed

python misalignment_check.py

## 3.4. graph features

-use make_graph_features_for_myproteins.py to extract graph features

python make_graph_features_for_myproteins.py

# 4. Run the model

For running the model, we need to the following steps

## 4.1. modifying DPFunc_pred.py
** DPFun_pred is already in DPFunc files, you need to find it and replace the line below:

checkpoint = torch.load('./save_models/{0}_{1}_{2}of{3}model.pt'.format(pre_name, ont, i_t_min, 3))

with this line

checkpoint = torch.load(
    './save_models/{0}_{1}_{2}of{3}model.pt'.format(pre_name, ont, i_t_min, 3), map_location='cpu')

## 4.2. go.obo

-Download go.obo and add it to ~/DPFunc/data directory

wget http://purl.obolibrary.org/obo/go.obo -O data/go.obo

-we also need to copy a file related to mf module of the model for bp and cc so that we can run all mf, bp, and cc predictions:

cp -f ./processed_file/graph_features/mf_test_whole_pdb_part0.pkl ./processed_file/graph_features/bp_test_whole_pdb_part0.pkl
cp -f ./processed_file/graph_features/mf_test_whole_pdb_part0.pkl ./processed_file/graph_features/cc_test_whole_pdb_part0.pkl

## 4.3. edit configure files to link them to our input files 

-within ~/DPFunc/configure directory, there are three .yaml files: mf.yaml, bp.yaml, and cc.yaml


-For all mf.yaml, bp.yaml, and cc.yaml file, replace the items below with their related lines

base:
  interpro_whole: ./processed_file/interpro_22369/{}.pkl

test:
  name: mytest
  pid_list_file: ./data/test_pid_list.pkl
  pid_go_file: ./processed_file/placeholder_go.txt
  pid_pdb_file: ./processed_file/graph_features/cc_test_whole_pdb_part0.pkl
  interpro_file: ./processed_file/cc_mytest_interpro_22369.pkl

## 4.4. Run the DPFunc predicaiton

To do the DPFunc prediciton, use Run_DPFunc.sb file

sbatch Run_DPFunc 

** After successfully completing the job, 12 .pkl files should be generated in the ~/DPFunc/results directory

# 5. Visualizing the outcome
DPFunc_model_mf_final.pkl, DPFunc_model_bp_final.pkl, and DPFunc_model_cc_final.pkl in the ~/DPFunc/results directory contains the resutls. Run the codes below to print the results in a human readable format:

run mf_print.py
run bp_print.py
run cc_print.py




# Contact

-If you have any issue/question please feel free to reach out to me at aminoroa@msu.edu.# DPFun_User_Guide
# DPFun_User_Guide
# DPFun_User_Guide
# DPFun_User_Guide
