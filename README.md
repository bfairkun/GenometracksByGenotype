# GenometracksByGenotype
Helper script to plot molQTLs. Normalize and group bigwigs by genotype.

## Dependencies
The python script requires a bunch of libraries:
- pyVCF
- pyBigWig
- numpy

Additionally, the intended use of the python script is to prepare bigwigs and track files for plotting with pyGenomeTracks.

All of the dependencies, and pyGenomeTracks, can be installed with the conda evnironment included:

```
conda create
```

## Usage

Use the main script with the help flag for more info

```
python
```

Also, included in `test_data/` is data from Grubert et al. Specifically, it is H3K4me3 ChIP-Seq data on chr19, and we will use it to plot a H3K4me3 QTL.

```
# Make a dir to store results
mkdir test_results

# Run 
python

# Run pyGenomeTracks

```

The following image is produced...
