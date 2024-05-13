# Installation
`pip install git+https://github.com/kburmak/plastic_mass_spectra.git`

# Library structure

## Database creation
`create_database` function creates `plastic_spectrum` database with three tables:
- `plastic` table, containing information about plastics' names;
- `spectra` table, containing information about spectra: plastic name, name of annotated fragment ion, \emph{m/z} ratio, peak intensity, (whether this ion generated on decoy step or not);
- `residue` table, containing information about fragment ions: plastic name, name of fragment ion, its \emph{m/z} ratio, and decoy status (whether this ion generated on decoy step or not).

## Generation of fragment ions
### Description
`generated_spectra` class contains three functions: `pss_gen` and `pfas` gen, which generate fragment ions for PSS and PFAS plastics respectively based on their chemical structure; and `plastic_to_database` function, which stores generated fragments in the `residue` database.

### Parameters
- `n_max`: maximum degree of polymerization.
- `n_mim`: minimum degree of polymerization.
- `max_so3na`: number of SO<sup><sub>3<\sub><\sup>Na losses.
- `max_so3`: number of SO<sup><sub>3<\sub><\sup> losses.

- `n_c`: number of carbons.
- `n_f`: number of fluorides.
- `n_k`: number of potassiums.
- `n_s`: number of sulphurs.
- `n_o`: number of oxygens.
- `n_h`: number of hydrogens.

### Methods
- `__init__`: initializes the class.

## Downloading experimental spectra, their annotation and plotting
### Description
`plastic_spectra` class contains five functions:
- `read_spectra` function reads spectra from `.ms1` files based on the retention time, and cleans it based on the threshold.
- `annotate_spectra` function annotates spectra based on the tolerance gap.
- `save_spectra` function saves spectra (annotated or not) with the proposed name.
- `load_spectra` function loads chosen spectra from the database.
- `plot_spectra` function plots spectra in three ways: experimental spectra, annotated experimental spectra with any random matches (decoy or not), and experimental spectra without random matches.

### Parameters
- `filename`: name of the file.
- `retention_time`: retention time of the spectra.
- `threshold`: threshold for the `find_peaks` function.
- `tolerance`: tolerance for matching theoretical and experimental m/z ratios.
- `name`: name of spectra to save in database and load from it.
- `title`: title for plots.

### Methods
- `__init__`: initializes the class.
- `__residue_search`: searches in the database based on m/z ratio and tolerance gap.
- `__clean`: describe status of match in the annotation (whether random match or not).
