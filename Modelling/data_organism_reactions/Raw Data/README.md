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

Raw Data contains all the raw unprocessed files

This folder should contain four files:

* all_genomes.txt
* ModelSEED-compounds-db.csv
* ModelSEED-reactions-db.csv
* table_with_reactions.txt

folder: locations
	=> contains raw files for locations


