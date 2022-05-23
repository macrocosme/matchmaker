MATCHMAKER
----------

MATCHMAKER is a python package to perform cross-matching of catalogs and parameter comparisons. 

USAGE
-----

Firstly, create a `Catalog` instance within `matchmaker.matchmaker.datasets` for each of your datasets. Example instances are provided.
The instance requires a `load_data` function that loads the dataset into a pandas Dataframe. Current APIs are provided for input files 
in `csv` and `fits` format. At least two datasets are required to perform cross-matching (or compare a unique dataset to itself).
For each dataset, the Catalog class requires, at minimum, how to access the columns corresponding to right ascension (RA) and declination (DEC).

Secondly, cross-matching is performed using `matchmaker.matchmaker.crossmatch.crossmatch`, specifying two `Catalog` inputs, 
spatial separation threshold for matching with corresponding physical unit, and wether or not to filter within an ellipse 
around the source centroid. Catalog masks can be provided to work on partial datasets.

Thirdly, evaluate chance allignment probabilities with one of `matchmaker.matchmaker.stats.chance_association` (circle) or 
`matchmaker.matchmaker.stats.chance_association_ab` (ellipse), specifying the source catalog `Catalog`, the name of the target `Catalog`, 
and internal indexes of matched sources, number of sources in target catalog and the area they cover (to get the target catalog source density). 

Results are stored under `Catalog(name='self').matches[Catalog(name='other').name]` and include a variety of information about the cross-section
(`mask` for source, matched indexes for target (`idx`), shortcuts (`mask_idx`=`numpy.where(mask)[0]`, `filtered_idx`=`idx[mask]`), `chance` association).

Loaded data in a Catalog object can be accessed as a `Pandas.Dataframe` instance, letting use all of pandas functionalities (`.loc`, `.iloc`, `.iterrows()`, `.isna()`). 
Further access parameters can be found via `Catalog().cols` to retreive dataframe column names programmatically (e.g. `Catalog().cols.<column_name>.label`). 
Similarly, columns units can be found under `Catalog().cols.<column_name>.unit`. 

