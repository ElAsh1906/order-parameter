# Lipid Order Parameter Calculation

This project contains a script for calculating the order parameters of lipid tails from molecular dynamics (MD) simulation data. The order parameter is a measure of the orientation of lipid tails relative to a reference direction, often the bilayer normal in lipid membranes. It provides insights into the structural properties and dynamics of lipid bilayers. The main features of this script are the ability to select lipids closest to the drug molecule and to process lipids with tails of different saturation.

## Description

The lipid order parameter $S_{CD}$ is calculated using the following equation:

$$S_{CD} = \left< \frac{3 \cos^2(\theta) - 1}{2} \right>$$

where $\theta$ is the angle between the lipid tail vector and the reference direction (in this case the bilayer normal), and $\left< \cdot \right>$ denotes an ensemble average.

For more details on lipid order parameters, see [this review article](https://doi.org/10.1021/acs.jctc.7b00643).

## Installation and usage

To run this script, you need to have Python 3 and [MDAnalysis](https://www.mdanalysis.org/) installed. You can install MDAnalysis using pip:

```sh
pip install MDAnalysis
```
Clone this repository to your local machine:

```sh
git clone https://github.com/ElAsh1906/order-parameter.git
cd order-parameter
```

#### Usage: 
`
./order_params.py [-h] [-s <.tpr>] [-f <.xtc>] [-o <.dat/.txt/...>] -l {popc,dopc,dmpc} [-d <distance>] [--nearest] [--refmol <resname>] [-bf <first frame>] [-b <start time>] [-ef <last frame>] [-e <last time>] [-step <step>] [-ts <time step>]
`
#### Options to specify input files are:

`-s` &emsp;&emsp; <.tpr> &emsp;&emsp; \(topol.tpr\)  

&emsp;&emsp;&emsp;&emsp;Input topology file.  

`-f` &emsp;&emsp; <.xtc> &emsp;&emsp; \(traj_comp.xtc\)  

&emsp;&emsp;&emsp;&emsp;Input trajectory file.  

    
#### Options to specify output files are:

`-o` &emsp;&emsp; <.dat/.txt/...> &emsp;&emsp; \(order_params.txt\)  

&emsp;&emsp;&emsp;&emsp;Output file for writing calculated order parameters.   


#### Other options are:

`-l, --lipids`  

&emsp;&emsp;&emsp;&emsp;Lipid model used in the simulation (choices: popc, dopc, dmpc) (required).   

`-d, --distance` &emsp;&emsp; \(5.0 Å\)  

&emsp;&emsp;&emsp;&emsp;Threshold distance in Angstroms for identifying nearest lipids.  

`--nearest`  

&emsp;&emsp;&emsp;&emsp;Enable calculation of nearest lipids.  

`--refmol` &emsp;&emsp; \(MOL\)  

&emsp;&emsp;&emsp;&emsp;Residue name of the reference molecule around which the environment is studied (default: MOL).  

`-bf` &emsp;&emsp; <first frame> &emsp;&emsp; \(0\)  

&emsp;&emsp;&emsp;&emsp;Frame number to start the analysis from.  

`-b` &emsp;&emsp;&ensp; <start time> &emsp;&emsp; \(0 ps\)  

&emsp;&emsp;&emsp;&emsp;Start time in picoseconds (ps) to begin the analysis.  

`-ef` &emsp;&emsp; <last frame> &emsp;&emsp; \(-1\)  

&emsp;&emsp;&emsp;&emsp;Frame number to end the analysis at.  

`-e` &emsp;&emsp;&ensp; <last time> &emsp;&emsp; \(-1\)  

&emsp;&emsp;&emsp;&emsp;End time in picoseconds (ps) to stop the analysis.  

`-step` &ensp; <step> &emsp;&emsp; \(1\)  

&emsp;&emsp;&emsp;&emsp;Interval between frames to analyze.  

`-ts` &emsp;&emsp; <time step> &emsp;&emsp;  

&emsp;&emsp;&emsp;&emsp;Time step interval in picoseconds (ps) for the analysis.  


## Output

Output file is text file with two columns separated by comma. Below is an example of the output format:

```
Saturated Unsaturated
0.085,0.117
0.099,0.146
0.175,0.173
...
```
Each value corresponds to the order parameter for a particular carbon in the chain, numbering starts from the first carbon atom after the carboxyl group.  

If your lipids have the same tails (DOPC, for example), output will be divided into two columns as well. 

## License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Contact

If you have any questions or need further assistance, feel free to contact elenayakush43@gmail.com
