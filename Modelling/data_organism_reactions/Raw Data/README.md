## The Raw data

### 0. Overview

File | Description
------------| -----------
all_genomes.txt | a file with the genome id and the full genome names
ModelSEED-compounds-db.csv | compound id and the full name of the compounds
ModelSEED-reactions-db.csv | reaction id and the full description of the reactions
table_with_reactions.txt | genome id and the reaction id associated 

### 1. Metabolic Information

The data was downloaded via the [Network-Based SEED API](http://blog.theseed.org/servers/). 

#### Installing RastShell and accessing serverscripts
* Installation Information can be found [here] (http://blog.theseed.org/servers/)
* Information on using the Server Scripts can be found [here](http://pubseed.theseed.org/sapling/server.cgi?pod=ServerScripts) 

#### Downloading Raw data

* Via **RastShell** all reactions per model (this might take a while...): 
```perl
svr_all_models|svr_reactions_in_model> table_with_reactions.txt
```  
* Via **Rastshell** all genome names and genome IDs : 
```perl
svr_all_genomes -complete > all_genomes.txt
```  
* All the reactions in the database: [here](http://seed-viewer.theseed.org/ModelSEEDdownload.cgi?biochemistry=1).

* All the compounds in the database: [here](http://seed-viewer.theseed.org/ModelSEEDdownload.cgi?biochemCompounds=1).

Note: If any of the links are broken, try: [http://seed-viewer.theseed.org/seedviewer.cgi?page=ModelView](http://seed-viewer.theseed.org/seedviewer.cgi?page=ModelView).


### 2. Co-occurence in ecosystems 

Information on the [locations](./locations) 
can be retrieved via this [link](http://mblnx-kallisto.uzh.ch:8888/microbial_coexistence/) (Chaffron et al.). 

Two files to download: 
- gg_otus_allseqs.tgz
- gg_allseqs2envo.tsv

## Bibliography

* Chaffron, Samuel, et al. "A global network of coexisting microbes from environmental and whole-genome sequence data." Genome research 20.7 (2010): 947-959.

