# MultiEM: An OpenMM-based software for the whole genome 3D structure reconstruction

## Installation
The anaconda environment can be found in `biosym.yml` file. The software mainly needs OpenMM o run.

## Usage
All the parameters of the model are saved within a `config.ini` file. Having specified the paramameters that you would like to use, you can run on terminal the following command,

```
python MultiEM.py -c config.ini
```

The sotware will return you a folder with the resulting structure, and some plots that show how compartments are distributed.
