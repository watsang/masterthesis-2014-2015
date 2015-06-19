Data retrieval
--------------
Finding the data is described in retrieval_data.txt.

Retrieving the data for constructing similarity metrics can be performed via the RAST-Server.

### Installing RastShell and accessing serverscripts
* Information can be found [here](http://blog.theseed.org/servers/). 
* Access the server scripts [here](http://pubseed.theseed.org/sapling/server.cgi?pod=ServerScripts). 

### Raw data

1. All reactions per model via **RastShell**  : 
```perl
svr_all_models|svr_reactions_in_model>all_reactions_in_table.txt
```  
2. All genome names and genome IDs via **Rastshell**: 
```perl
svr_all_genomes -complete > all_genomes.txt
```  
3. All the reactions in the database: 
Via this [link](seed-viewer.theseed.org/ModelSEEDdownload.cgi?biochemistry=1).

4. All the compounds in the database: 
Via this [link](seed-viewer.theseed.org/ModelSEEDdownload.cgi?biochemCompounds=1).

Note: If any of the links are broken, try: [http://seed-viewer.theseed.org/seedviewer.cgi?page=ModelView](http://seed-viewer.theseed.org/seedviewer.cgi?page=ModelView).

Data preprocessing
-------------------
The data preprocessing files are in the directory data_preprocessing. 
These files are unnecessary for the ensuing analysis. 

The raw data for the organism reactions are in the directory data_organism_reactions.

Data analysis
------------------
An initial pca analysis was performed in the directory pca. 
For this, you will need the data in the folder data_organism_reactions.