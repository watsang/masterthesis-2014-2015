### Experimental Design
The R-file for generating the design can be found in [Re_preparation co-culture](./Re_ preparation co-culture). 

The raw output of the design is design.csv. The R-file was changed to a computational lighter version which generates only the relevant variables. 

#### Aim of the design

The research questions that were tackled in my thesis mainly deal with synthetic microbial ecology. As pathogen increasingly gain antibiotic resistance, it becomes more interesting to look at the underlying factors of community structure to manipulate a community's resistance to invasion. One key factor that is explored is community evenness, i.e. the relative abundance of species in a community. It is reported in literature that communities with low evenness -- meaning strong domination by one species -- is more fragile than communities with high evenness -- i.e. communities where most species are equally abundant. 

#### Formulating a null-hypothesis

A hypothesis was formulated that communities with high evenness are more resistant to invasion contrary to communities with low evenness. Lab experiments were set up where 194 synthetic communities consisting of 10 bacteria were assembled on a varying range of evenness. To find the optimal design 1000000 communities were simulated, after which a 94 were randomly selected and 100 according to a certain stratification. 

#### Description variables of design.
```b1``` to ```b10```: The total cellconcentration of the ten bacteria in the synthetic community
```evenness```: The evenness calculated with pielou's evenness based on the cellconcentration of the ten bacteria ```b1``` to ```b10```. 

