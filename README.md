# drrglobal
DRRGlobal: recovering weak phases from global seismograms 

## Reference
If you find this package useful, please do not forget to cite the following paper.

    Chen and Chen, 2024. DRRGlobal: uncovering the weak phases from global seismograms using the damped rank-reduction method, under review.
    
BibTeX:
	
	@article{drrglobal,
	  author={Wei Chen and Yangkang Chen},
	  title = {Uncovering the weak phases from global seismograms using the damped rank-reduction method},
	  journal={TBD},
	  year=2023,
	  volume=X,
	  issue=X,
	  pages={under review},
	  doi={XXX},
	}

-----------
## Copyright
    Authors of the drrglobal paper, 2021-present
-----------

## License
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)   

-----------

## Install
Using the latest version

    git clone https://github.com/chenyk1990/drrglobal
    cd drrglobal
    addpath(genpath('./')); #in Matlab command line
    
-----------
## Examples
    The "main" directory contains all runable scripts test_figNO.m to reproduce all figures in the paper.
 
-----------
## Dependence Packages
* MATdrr (https://github.com/chenyk1990/MATdrr)
    
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

4. The directory [subroutines](https://github.com/chenyk1990/drrglobal/tree/main/subroutines) stores all the required subroutines. 

5. The current version is based on Matlab. Future versions may also support Python and be optimized regarding computational efficiency. 

6. All figures (except for fig1, which is a schematic plot, and fig11, which is based on Madagascar and powerpoint) in the drrglobal paper are in the following directory for a quick look (https://github.com/chenyk1990/drrglobal/tree/main/gallery/). 

7. The label fonts are different across different Matlab versions and platforms (Linux, Max). The figures presented in the paper are from Mac-Pro Matlab 2022b. 

-----------
## Gallery
The gallery figures of the drrglobal package can be found at
    https://github.com/chenyk1990/gallery/tree/main/drrglobal
Each figure in the gallery directory corresponds to a test_figNO.m script in the "main" directory with the exactly the same file name (figNO.png).

These gallery figures are also presented below. 

Following Figures Generated by [test_real.m](https://github.com/chenyk1990/drrglobal/tree/main/test_real.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/drrglobal/real_original.png' alt='Slicing' width=960/>
<img src='https://github.com/chenyk1990/gallery/blob/main/drrglobal/real_rr.png' alt='Slicing' width=960/>
