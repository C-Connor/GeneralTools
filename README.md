# GeneralTools
Collection of custom scripts for performing various tasks.

## FeaturesFromGBK.py
Python script to extract features from GBK file. Type of feature can be viewed using the `--list_features` option.

## CDSToAminoAcid.py
Python script to translate CDS multifasta file to peptide multifasta file.

## GeneListToMultifasta.py
Extracts sequences of interest from multifasta file. Takes a list of sequence identifiers, one per line, and extracts corresponding sequences from the multifasta file.

## ExtractGeneClustersToIndividualFasta.py
Extracts the Creates individual multifasta files for each cluster in the given Roary/CD-HIT cluster file. Requires the individual gene sequences for each of the isolates in the cluster file e.g. output from prokka.