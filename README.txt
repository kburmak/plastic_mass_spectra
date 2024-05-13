# Installation
`pip install git+https://github.com/kburmak/plastic_mass_spectra.git`

# Library structure

## Database creation
`create_database` function creates `plastic_spectrum` database with three tables:
1. \emph{plastic} table, containing information about plastics' names;
2. \emph{spectra} table, containing information about spectra: plastic name, name of annotated fragment ion, \emph{m/z} ratio, peak intensity, (whether this ion generated on decoy step or not);
3. \emph{residue} table, containing information about fragment ions: plastic name, name of fragment ion, its \emph{m/z} ratio, and decoy status (whether this ion generated on decoy step or not).

## Generation of fragment ions
### Description
This class contains two functions: `pss_gen` and `pfas` gen, which generate 
