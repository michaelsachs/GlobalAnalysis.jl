
## File format

Experimental data should be formatted as a .csv file, where:
- The first column (second entry onwards) is the spectral dimension, e.g. wavelength/wavenumbers/energy/etc.
- The first row (second entry onwards) is the temporal dimension
- The first value of the dataset (first column, first row) is unused; for example, it can be set to `0`
- The rest of the file are amplitudes of the 2-dimensional spectro-temporal data

An example dataset can be found in `\data\testData_first_order_seq.csv`.

## Data import

Data is read from .csv files using `importData(path)`, where `path` can point either to a single .csv file or a folder.  
- If `path` points to a folder, all .csv files contained in it are imported
- If `path` points to a .csv file, only this file is imported

In both cases, `importData` returns a [StructArray](https://juliaarrays.github.io/StructArrays.jl/stable/), containing one or more elements depending on how many .csv files were imported, which can be accessed as follows:

```julia
# define .csv file containing experimental data
file = raw"C:\GlobalAnalysis.jl\data\testData_mixed_order_par.csv"
# import data, pick first dataset
data = importData(file)[1]

# spectral dimension
s = data.x
# time dimension
t = data.y
# 2D data
d = data.z
```