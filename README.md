# sigrecover
SigRecover: Recovering signal from noise in distributed acoustic sensing data processing

## Reference
If you find this package useful, please do not forget to cite the following paper.

    Chen, Y., 2023. SigRecover: Recovering signal from noise in distributed acoustic sensing data processing, under review.
    
BibTeX:
	
	@article{sigrecover,
	  author={Yangkang Chen},
	  title = {SigRecover: Recovering signal from noise in distributed acoustic sensing data processing},
	  journal={TBD},
	  year=2023,
	  volume=X,
	  issue=X,
	  pages={under review},
	  doi={XXX},
	}

-----------
## Copyright
    Authors of the sigrecover paper, 2021-present
-----------

## License
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)   

-----------

## Install
Using the latest version

    git clone https://github.com/chenyk1990/sigrecover
    cd sigrecover
    addpath(genpath('./')); #in Matlab command line
    
-----------
## Examples
    The "main" directory contains all runable scripts test_figNO.m to reproduce all figures in the paper.
 
-----------
## Dependence Packages
* MATseisdl (https://github.com/chenyk1990/MATseisdl)
    
-----------
## Development
    The development team welcomes voluntary contributions from any open-source enthusiast. 
    If you want to make contribution to this project, feel free to contact the development team. 

-----------
## Contact
    Regarding any questions, bugs, developments, collaborations, please contact  
    Yangkang Chen
    chenyk2016@gmail.com

-----------
## NOTES:
 
1. To run the reproducible scripts (test_xxxx.m), please first download the required package from: https://github.com/chenyk1990/MATseisdl. 

2. Please put the MATseisdl package in the main directory and add its sub-directories to the Matlab toolbox path. 

3. The scripts beginning with "test_" are runnable scripts.

4. The directory [subroutines](https://github.com/chenyk1990/sigrecover/tree/main/subroutines) stores all the required subroutines. 

5. The current version is based on Matlab. Future versions may also support Python and be optimized regarding computational efficiency. 

6. All figures (except for fig1, which is a schematic plot, and fig11, which is based on Madagascar and powerpoint) in the sigrecover paper are in the following directory for a quick look (https://github.com/chenyk1990/sigrecover/tree/main/gallery/). 

7. The label fonts are different across different Matlab versions and platforms (Linux, Max). The figures presented in the paper are from Mac-Pro Matlab 2022b. 

-----------
## Gallery
The gallery figures of the sigrecover package can be found at
    https://github.com/chenyk1990/gallery/tree/main/sigrecover
Each figure in the gallery directory corresponds to a test_figNO.m script in the "main" directory with the exactly the same file name (figNO.png).

These gallery figures are also presented below. 

Figure 2 Generated by [test_jdas.m](https://github.com/chenyk1990/sigrecover/tree/main/test_jdas.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/sigrecover/jdas.png' alt='Slicing' width=960/>
<img src='https://github.com/chenyk1990/gallery/blob/main/sigrecover/jdas_z.png' alt='Slicing' width=960/>
