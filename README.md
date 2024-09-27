# eos_workflow
A workflow for eos calculations

## Installation Prerequisites

### ABINIT
The installation of ABINIT can be referred to their official website

### MongoDB
Details can be referred to the website of atomate2

### Configure files for atomate2 and jobflow

#### i. configure setting files
atomate2.yaml
```
# ABINIT
ABINIT_MPIRUN_CMD: "mpirun"
ABINIT_CMD: PATH_TO/src/98_main/abinit
ABINIT_MRGDDB_CMD: PATH_TO/src/98_main/mrgddb
ABINIT_ANADDB_CMD: PATH_TO/src/98_main/anaddb
ABINIT_ABIPY_MANAGER_FILE: $HOME/.abinit/abipy/manager.yml
```

jobflow.yaml
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

manager.yml, abipy configuration file
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
pip install -e .
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
