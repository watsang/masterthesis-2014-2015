### Experimental Design
The R-file for generating the design can be found in [Re_preparation co-culture](./Re_ preparation co-culture). 

The raw output of the design is design.csv. The R-file was changed to a computational lighter version which generates only the relevant variables. 

#### Aim of the design

The research questions that were tackled in my thesis mainly deal with synthetic microbial ecology. As pathogen increasingly gain antibiotic resistance, it becomes more interesting to look at the underlying factors of community structure to manipulate a community's resistance to invasion. One key factor that is explored is community evenness, i.e. the relative abundance of species in a community. It is reported in literature that communities with low evenness -- meaning strong domination by one species -- is more fragile than communities with high evenness -- i.e. communities where most species are equally abundant. 

#### Formulating a null-hypothesis

A hypothesis was formulated that communities with high evenness are more resistant to invasion contrary to communities with low evenness. Lab experiments were set up where 194 synthetic communities consisting of 10 bacteria were assembled on a varying range of evenness. To find the optimal design 1000000 communities were simulated, after which a 94 were randomly selected and 100 according to a certain stratification. 

#### Description variables of design.
* ```b1``` to ```b10```: The total cells of the ten bacteria in the synthetic community. Values between 10^3 to 10^7 cells per volume cryovial. 
* ```evenness```: The evenness calculated with pielou's evenness based on the cells of the ten bacteria ```b1``` to ```b10```. 
* ```celaantal```: the total cell count in the specific community. The sum of all the cells for the ten bacteria. 
* ```gini_thas```:  the [gini-coefficient](https://en.wikipedia.org/wiki/Gini_coefficient), a measure of evenness often used in economy and ecology. This gini is calculated with the formula in [Environmental conditions and community evenness determine the outcome of biological invasion](http://www.nature.com/ncomms/journal/v4/n1/full/ncomms2392.html).
* ```gini_wittebolle```: the [gini-coefficient](https://en.wikipedia.org/wiki/Gini_coefficient) calculated with the formula in [Initial community evenness favours functionality under selective stress](users.ugent.be/~wverstra/research/Labmet%20publicatie%20922%20Wittebolle%20...%20Nature%20458,%20623-626.pdf)
* ```v1``` to ```v10```: The total volume of the corresponding ten bacteria in the synthetic community to be pipetted from the resp. concentration values ```c1``` to ```c10```. These variables have no meaning whatsoever for the analysis and are only created to make working tables for conducting the lab experiments. 
* ```c1``` to ```c10```: The concentration corresponding with the ten bacteria from which to pipet.These variables have no meaning whatsoever for the analysis and are only created to make working tables for conducting the lab experiments.


