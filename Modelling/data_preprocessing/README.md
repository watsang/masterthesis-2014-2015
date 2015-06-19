## The Raw data

The data was downloaded via the [Network-Based SEED API](http://blog.theseed.org/servers/). 

### Installing RastShell and accessing serverscripts
* Installation Information can be found [here] (http://blog.theseed.org/servers/)
* Information on using the Server Scripts can be found [here](http://pubseed.theseed.org/sapling/server.cgi?pod=ServerScripts) 

### Downloading Raw data

* Via **RastShell** all reactions per model (this might take a while...): 
```perl
svr_all_models|svr_reactions_in_model>all_reactions_in_table.txt
```  
* Via **Rastshell** all genome names and genome IDs : 
```perl
svr_all_genomes -complete > all_genomes.txt
```  
* All the reactions in the database: [here](seed-viewer.theseed.org/ModelSEEDdownload.cgi?biochemistry=1).

* All the compounds in the database: [here](seed-viewer.theseed.org/ModelSEEDdownload.cgi?biochemCompounds=1).

* The information on the locations: [here](http://mblnx-kallisto.uzh.ch:8888/microbial_coexistence/)

Note: If the links for the reactions and compounds are broken, try: [http://seed-viewer.theseed.org/seedviewer.cgi?page=ModelView](http://seed-viewer.theseed.org/seedviewer.cgi?page=ModelView).

## The Tidy data
The data preprocessing files are in the directory data_preprocessing. 
These files are unnecessary for the ensuing analysis. 

The raw data for the organism reactions are in the directory data_organism_reactions.

The data was preprocessed in the file:

split_into_reactions.ipynb