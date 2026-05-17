# FVT2DTHERMAL
FVT2DTHERMAL is an open-source MATLAB® implementation for verifying the temperature field in a two-dimensional homogeneous material with thermal conductivity 
𝑘, subject to imposed temperatures on its four sides (left, bottom, right, and top) using the Finite-Volume Theory (FVT).

# Table 1: Inputs parameters declaration
| Name | Description | Value |
|-------------|-------------|-------------|
| H | Height of the solid | 1,2,3,... |
| L | Width of the solid | 1,2,3,... |
| nx1 | Number of subvolumes along the x1 direction | 30, 40, 50, ... |
| nx2 | Number of subvolumes along the x2 direction | 30, 40, 50, ... |
| k | Thermal conductivity of the material | User-defined |
| fx1 | Subvolume selected along the x1 direction for plotting the temperature field along the x2 axis | Example |
| fx2 | Subvolume selected along the x2 direction for plotting the temperature field along the x1 axis | Example |
| tempFB | Temperature of the bottom faces | User-defined |
| tempFR | Temperature of the right faces | User-defined |
| tempFT | Temperature of the top faces | User-defined |
| tempFL | Temperature of the left faces | User-defined |

# Example: px1 and px2:
To exemplify the range of values for fx1 and fx2, we will use a discretized mesh of 50 x 50 subvolumes. In this case, the faces for fx1 are numbered in the following sequence: 1, 2, 3, 4, ..., 50, while the faces for fx2 follow the progression: 1, 51, 101, 151, ..., 2501.

# AUTHORS' CONTACT
Profº. Dr. Profº Dr. Márcio André Araújo Cavalcante - marcio.cavalcante@ceca.ufal.br  
Profº. Dr. Romildo dos Santos Escarpini Filho - romildo.escarpini@penedo.ufal.br  
PhD Student Diogo Tiago dos Santos - diogo.santos@ctec.ufal.br

