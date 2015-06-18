#### Statistical analysis

A hypothesis was formulated that communities with high evenness are more resistant to invasion contrary to communities with low evenness. The lab experiments were conducted and in each synthetic community per experiment one pathogen was introduced. 

Four pathogen are selected namely LMG7866, LMG3203, LMG7878 and LMG2954. The pathogen information can be found back in the file [human pathogens III def.xlsx](../pathogenlist). 

Let's go ahead and test this!

#### Description variables of design.

##### The variables in the columns
* ```b1``` to ```b10```: The total cells of the ten bacteria in the synthetic community. Values between 10^3 to 10^7 cells per volume cryovial. 
* ```evenness```: The evenness calculated with pielou's evenness based on the cells of the ten bacteria ```b1``` to ```b10```. 
* ```celaantal```: the total cell count in the specific community. The sum of all the cells for the ten bacteria. 
* ```gini_thas```:  the [gini-coefficient](https://en.wikipedia.org/wiki/Gini_coefficient), a measure of evenness often used in economy and ecology. This gini is calculated with the formula in [Environmental conditions and community evenness determine the outcome of biological invasion](http://www.nature.com/ncomms/journal/v4/n1/full/ncomms2392.html).
* ```gini_wittebolle```: the [gini-coefficient](https://en.wikipedia.org/wiki/Gini_coefficient) calculated with the formula in [Initial community evenness favours functionality under selective stress](users.ugent.be/~wverstra/research/Labmet%20publicatie%20922%20Wittebolle%20...%20Nature%20458,%20623-626.pdf)
* ```v1``` to ```v10```: The total volume of the corresponding ten bacteria in the synthetic community to be pipetted from the resp. concentration values ```c1``` to ```c10```. These variables have no meaning whatsoever for the analysis and are only created to make working tables for conducting the lab experiments. 
* ```c1``` to ```c10```: The concentration corresponding with the ten bacteria from which to pipet.These variables have no meaning whatsoever for the analysis and are only created to make working tables for conducting the lab experiments.
