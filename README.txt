# RNA_assessment_V2

This git includes:

- RNA structure normalization tool used in RNA-Puzzles.
- the RNA 3D structure comparison metrics used in RNA-Puzzles assessment, including RMSD (all atom), P value, Interaction Network Fidelity and Deformation Index.
- a Dockerfile for creating a Docker image in which the tools present in this git can be used immediately without having to install the dependencies on your personal machine
- the rmsd_di.py program for calculating RMSD and Deformation Index (DI) values ​​using RNA PDB files

## Dockerfile

To generate the Docker image from the provided Dockerfile, use the following command line when you are on this repository:
```
docker build -t rnaa_v2 .
```

From this image, you can now run your container using:
```
docker run -it --rm rnaa_v2 /bin/bash
```

**Warning :** the `--rm` option causes the container and its contents to be deleted when it is stopped

Once the container is running, a volume is created at the pathway `/var/lib/docker/volumes`. It allows the communication of files and directories between the host machine and the containe.
For example, to transfer a folder from the host machine to the container, you will need to enter in a second terminal a command like:
```
sudo cp -r /pathway/of/my/folder/ /var/lib/docker/volumes/random_volume_name/_data/folder
```

*Help :* `random_volume_name` *can be find using* `sudo ls /var/lib/docker/volumes`.

Your folder will be copied in the directory named `data_DI` of your container.
You can also recover the files or folders contained in this container by copying them to the `data_DI` directory, then executing the following type of command in the second terminal :
```
sudo cp -r /var/lib/docker/volumes/random_volume_name/_data/folder /pathway/of/my/host/machine/folder/ 
```

For more information about the use of Docker :
- [How to install Docker](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-20-04-fr)
- [Start to use and understand Docker](https://takacsmark.com/dockerfile-tutorial-by-example-dockerfile-best-practices-2018/)

## Installation on your personal machine
If you don't want to use Docker, you can also install the package:
`git clone https://github.com/RNA-Puzzles/RNA_assessment.git`
`cd RNA_assessment`
`python setup.py install`

This package depends on `biopython`.
To install: `pip install biopython`

## rmsd_di.py
This program is used to calculate the RMSD and the Deformation Index between a native structure and its experimental structure, or a folder of experimental structures.

### Examples of use
**Native vs Experimental (1 vs 1)**
```
python3 rmsd_di.py -n data_DI/1Z43/1Z43.pdb -e data_DI/1Z43/EXP_1Z43.pdb
```

The `-n` option is followed by the PDB of the native structure (which is required), and the `-e` option is followed by the PDB of its experimental structure.

**Native vs Folder of experimental (1 vs many)**
Basically, you can simply run this program in the folder containing experimental PDB structures to have RMSD and DI of all files in this folder, always specifying the native structure to use:
```
python3 rmsd_di.py -n data_DI/1Z43/1Z43.pdb
```

If you are not in this folder, you will have to specify the pathway of this folder with the option `-d`:
```
python3 rmsd_di.py -n data_DI/1Z43/1Z43.pdb -d data_DI/1Z43/
```

**More options**
- During the execution of the program, files like "1Z43.pdb.mcout" are created and then deleted by default once used. It is however possible to keep these files by adding the function `--view_mcout` to your command line.
- In the same way, it is possible to display all the values ​​that the program can calculate by adding the option `--all_values` (values : RMSD, DI, INF_ALL, INF_WC, INF_NWC, INF_STACK).
