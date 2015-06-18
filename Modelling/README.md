Data retrieval
--------------
Finding the data is described in retrieval_data.txt.

Retrieving the data for constructing similarity metrics can be performed via the RAST-Server.

Via the RastShell retrieve all reactions per model: 
```perl
svr_all_modells|svr_reactions_in_model>all_reactions_in_table.txt
```

Data preprocessing
-------------------
The data preprocessing files are in the directory data_preprocessing. 
These files are unnecessary for the ensuing analysis. 

The raw data for the organism reactions are in the directory data_organism_reactions.

Data analysis
------------------
An initial pca analysis was performed in the directory pca. 
For this, you will need the data in the folder data_organism_reactions.