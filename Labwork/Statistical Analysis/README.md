#### Statistical analysis

A hypothesis was formulated that communities with high evenness are more resistant to invasion contrary to communities with low evenness. The lab experiments were conducted and in each synthetic community per experiment one pathogen was introduced. 

Four pathogen are selected namely LMG7866, LMG3203, LMG7878 and LMG2954. The pathogen information can be found back in the file [human pathogens III def.xlsx](../pathogenlist). Four each pathogen the colonies are counted after 48 hours of invasion in the synthetic community. 

#### Description variables of design.

##### The variables in the columns
* ```id```: indicates which of the 194 communities is invaded. 
* ```b1``` to ```b10```: The total cells of the ten bacteria in the synthetic community. Values between 10^3 to 10^7 cells per volume cryovial. 
* ```evenness```: The evenness calculated with pielou's evenness based on the cells of the ten bacteria ```b1``` to ```b10```. 
* ```cellcount```: the total cell count in the specific community. The sum of all the cells for the ten bacteria. 
* ```gini_thas```:  the [gini-coefficient](https://en.wikipedia.org/wiki/Gini_coefficient), a measure of evenness often used in economy and ecology. This gini is calculated with the formula in [Environmental conditions and community evenness determine the outcome of biological invasion](http://www.nature.com/ncomms/journal/v4/n1/full/ncomms2392.html).
* ```gini_wittebolle```: the [gini-coefficient](https://en.wikipedia.org/wiki/Gini_coefficient) calculated with the formula in [Initial community evenness favours functionality under selective stress](users.ugent.be/~wverstra/research/Labmet%20publicatie%20922%20Wittebolle%20...%20Nature%20458,%20623-626.pdf)
* ```LMG7866```: Colonies counted of Escherichia fergusonii after invasion and plating on LB agar. 
* ```LMG3203```: Colonies counted of Raoultella terrigena after invasion and plating on LB agar. 
* ```LMG7878```: Colonies counted of Proteus inconstans after invasion and plating on LB agar. 
* ```LMG2954```: Colonies counted of Proteus mirabilis after invasion and plating on LB agar. 
* ```replica```: indicates the replica number. Each of the 194 communities is replicated two times. 

#### The statistics
The data is analysed in the file comparison_freqvsbayes.R. In a first part a frequentist approach is used. In a second part of the file an attempt was made to apply Bayes. 
###### Poisson
First a Poisson regression is used. But the model for each of the four pathogen are shown to be overdispersed. 
###### Negative-Binomial Regression
A negative binomial regression model can account for the overdispersion
###### Checking assumptions
Homoscedasticity is checked.


Let's go ahead and test our hypothesis!