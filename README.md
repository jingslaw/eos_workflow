# eos_workflow
A workflow for eos calculations

## installation
#### clone the code
```
git clone https://github.com/jingslaw/eos_workflow.git
```
#### create a conda env *jobflow*
```
conda create --name jobflow python=3.10
conda activate jobflow
```
#### install
```
pip install .
```

## Basical settings

### MongoDB
Details can be referred to the website of atomate2

### Abipy
Details can be referred to the website of abipy

### Configure files for atomate2, jobflow, and abipy
```
conda activate jobflow
pip show atomate2
```
You will obtain the path to atomate2: PATH_TO_ATOMATE2, like "/home/wjing/miniconda3/envs/phonon/lib/python3.10/site-packages". \
Then
```
cd PATH_TO_ATOMATE2/atomate2
```
and create a folder to store configuration files
```
mkdir config
cd config
```
The next step is to setting three configure files.

#### i. configure setting files
1. manager.yml \
abipy configuration file, usually store at "$HOME/.abinit/abipy". More details please referr to their website \
a template is here
```
qadapters:
    # List of qadapters objects
    - priority: 1
      queue:
        qtype: shell
        qname: localhost
      job:
        mpi_runner: mpirun
        # source a script to setup the environment.
        #pre_run: "source ~/env.sh"
      limits:
        timelimit: 1:00:00
        max_cores: 2
      hardware:
         num_nodes: 1
         sockets_per_node: 1
         cores_per_socket: 2
         mem_per_node: 4 Gb
```

2. atomate2.yaml \
a template as follows.
```
# ABINIT
ABINIT_MPIRUN_CMD: "mpirun"
ABINIT_CMD: PATH_TO/src/98_main/abinit
ABINIT_MRGDDB_CMD: PATH_TO/src/98_main/mrgddb
ABINIT_ANADDB_CMD: PATH_TO/src/98_main/anaddb
ABINIT_ABIPY_MANAGER_FILE: $HOME/.abinit/abipy/manager.yml
```
Notice:\
ABINIT_MPIRUN_CMD is only required when you use mpi to run abinit, if not, just delete it.\
ABINIT_ABIPY_MANAGER_FILE show the manager.yml required by Abipy

3. jobflow.yaml \
Here is a template. \<\<DB_NAME\>\>, \<\<HOSTNAME\>\>, \<\<PORT\>\>, \<\<USERNAME\>\>, and \<\<PASSWORD\>\> should be replaced by your own MongoDB
```
JOB_STORE:
  docs_store:
    type: MongoStore
    database: <<DB_NAME>>
    host: <<HOSTNAME>>
    port: <<PORT>>
    username: <<USERNAME>>
    password: <<PASSWORD>>
    collection_name: outputs
  additional_stores:
    data:
      type: GridFSStore
      database: <<DB_NAME>>
      host: <<HOSTNAME>>
      port: <<PORT>>
      username: <<USERNAME>>
      password: <<PASSWORD>>
      collection_name: outputs_blobs
```

#### ii. export environment variables

open ~/.bashrc and add the following:
```
export ATOMATE2_CONFIG_FILE="PATH_TO/atomate2.yaml"
export JOBFLOW_CONFIG_FILE="PATH_TO/jobflow.yaml"
```

then:
```
source ~/.bashrc
```
## Install Pseudopotentials
We use abipy to install pseudopotentials.\
The following line will show some avialiable pseudopotential families.
```
abips.py avail
```
Then we choose "ONCVPSP-PBE-SR-PDv0.4" as an example to install it. 
```
abips.py install ONCVPSP-PBE-SR-PDv0.4
```
ONCVPSP means pseudopotentials are norm-conserving. 
PBE is the type of exchange-correlation functional. SR means the pesudopotentials are scalar-relativistic.
v0.4 is the version of this pseudopotentials family.


## An Example for O.psp8
An example is in eos_workflow/src/tests/ \
pseudo O.psp8 is from "ONCVPSP-PBE-SR-PDv0.4:standard", a standard
PseudoTable, which installed in abipy:
~/.abinit/pseudo/ONCVPSP-PBE-SR-PD/standard

### run locally
```
python run_locally.py
```
Many folders will create to store calculation files and final results.
You will find a file named "eos_fitting_results.json" at the last created folder. If everything is fine, the result is like "eos curve is bad". \
This is what we expect, because we use a very small ecut as testing the workflow here.  

### run the workflow on a remote cluster

There are two ways to realize it: \
a. run the workflow locally on the remote cluster.\
In this case, you need to write a bash script to correctly submit your job to the calculation node rather than your home node.\
In that bash script, you could still use something like "python run_locally.py"\

b. run the workflow locally on your own PC, and let jobflow remote submit it to the cluster\
In this case, we run workflows directly on our local PC.

#### i. install eos_workflow package on the remote cluster
same steps as above to install eos_workflow, set configure files on cluster
#### ii. install eos_workflow package on PC
#### iii. install jobflow_remote package on PC
When eos_workflow is installed on conda env jobflow, we
follow the introduction of jobflow_remote:
```
pip install jobflow-remote
```
Then, initial setup configuration:
```
jf project generate eos_workflow
```
In addition, file ~/.jfremote.yaml should be created. \
~/.jfremote.yaml with one line:
```
project: eos_workflow
```
Finally, create and configure the eos_workflow.yaml file
in the folder ~/.jfremote

If this file is correctly configured, type the following line
```
jf project list
```
and we will see a project named eos_workflow in green. If that "eos_workflow" is white, you need to modify your jobflow_remote settings.

#### iv. submit a remote flow
Here is an example in eos_workflow/src/tests/
```
python submit_remote.py
```
The submitted jobs can be referred by:
```
jf job list
```
