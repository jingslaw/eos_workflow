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
You will obtain the path to atomate2: PATH_TO_ATOMATE2, like "/home/wjing/miniconda3/envs/phonon/lib/python3.10/site-packages".
Then
```
cd PATH_TO_ATOMATE2/atomate2
```
and create a folder to store configuration files
```
mkdir config
cd config
vim atomate2.yaml
vim jobflow.yaml
```

#### i. configure setting files
manager.yml, abipy configuration file, usually store at $HOME/.abinit/abipy, more details please referr to their website
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

atomate2.yaml
a template as follows.
ABINIT_MPIRUN_CMD is only required when you use mpi to run abinit, if not, just delete it.
ABINIT_ABIPY_MANAGER_FILE show the manager.yml required by Abipy
```
# ABINIT
ABINIT_MPIRUN_CMD: "mpirun"
ABINIT_CMD: PATH_TO/src/98_main/abinit
ABINIT_MRGDDB_CMD: PATH_TO/src/98_main/mrgddb
ABINIT_ANADDB_CMD: PATH_TO/src/98_main/anaddb
ABINIT_ABIPY_MANAGER_FILE: $HOME/.abinit/abipy/manager.yml
```

jobflow.yaml
Here is a template. <<DB_NAME>>, <<HOSTNAME>>, <<PORT>>, <<USERNAME>>, and <<PASSWORD>> should be replaced by your own MongoDB
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

## An Example for O.psp8
An example is in eos_workflow/src/tests/ \
pseudo O.psp8 is from "ONCVPSP-PBE-SR-PDv0.4:standard", a standard
PseudoTable, which installed in abipy:
~/.abinit/pseudo/ONCVPSP-PBE-SR-PD/standard

### run locally
```
python run_locally.py
```

### run the workflow on a remote cluster

#### i. install eos_workflow package on the remote cluster
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

If this file is correctly configured, we have
```
python submit_remote.py
```
The submitted jobs can be referred by:
```
jf job list
```
