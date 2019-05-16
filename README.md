Repository for my projects during my secondment at Brandes University

### Structure

There is one central Snakefile in the root, which calls other Snakefiles under shared/code/snakefiles/.

Also, there is one config.yaml configuration file, which is currently configured to work with data at
the Rajewsky lab server, murphy.

Upon running snakemake one should use the ```--use-cores``` flag, with at least 16 cores, and also the ```--conda``` flag, as some 
scripts will use a conda environment with python2.

The script will generate a folder structure like ```<person>/<project>/data```. Every file generated will be found under this. As a result,
the ```data``` folders are always in the .gitignore'
