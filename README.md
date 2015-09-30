Genotype Phenotype Correlation Tool
========

We present GePhCorT (Genotype Phenotype Correlation Tool), an algorithm that relies on finding the optimal states of the characters over a phylogeny, and then using randomization testing to assess the significance of the correlation between the genotype and phenotype.

# Download

* Clone the repository. 
```
git clone https://github.com/Farhat/gephcort.git
```

* Enter the directory.

```
cd gephcort
```
# Dependencies

* Python 2.7
* R
* Redis

# Setup

* Install redis by typing the following commands into your terminal.
```
wget http://download.redis.io/redis-stable.tar.gz
tar xvzf redis-stable.tar.gz
cd redis-stable
make
make install
```
***

* Run the python setup file.
```
python setup.py install
```
***

* Start the redis server. 
```
redis-server
```
***

* Run the `resurrect.R` script. The generic command format is:
```
Rscript resurrect.R <sequence_file> <newick_tree_file> <fasta/phylip> <resurrect_output_file>
```
To try out the script on the example data, use the following command.
```
Rscript resurrect.R example/example.phy example/example.newick phylip resurrect_output
```
You can view the sample output by opening the `resurrect_output` file.

***

* Run the `reanimate.py` script. The generic command format is:
```
python reanimate.py -s <seq_file> -t <tree_file> -f <format(fasta/phylip)> -i <phen_iterations> -p <phen_file> -r <ressurect_output_file> -o <output_file> -c <cores>
```
To try out the script on the example data, use the following command.

```
python reanimate.py -s example/example.phy -t example/example.newick -f phylip -i 100 -p example/example_phenotype.txt -o gephcort_output -r example/example_resurrect.dat -c 1
```

You can view the sample output by opening the `gephcort_output` file.

