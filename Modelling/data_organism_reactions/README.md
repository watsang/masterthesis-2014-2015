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

Note: If any of the links are broken, try: [http://seed-viewer.theseed.org/seedviewer.cgi?page=ModelView](http://seed-viewer.theseed.org/seedviewer.cgi?page=ModelView).

## The Tidy data
The data preprocessing files are in the directory data_preprocessing. 
These files are unnecessary for the ensuing analysis. 

The raw data for the organism reactions are in the directory [Raw Data](./Raw Data).

The raw data of the reactions per organism can be found in directory [Raw Reaction Matrices](./Raw Reaction Matrices).

The [Stoichiometric Matrices](./Stoichiometric Matrices) folder contains all the processed Stoichiometric Matrices.
These are already processed. You can look up the reactions in csv with reactions.
The compounds can be replaced with the real compound name.
