# wildflower-traits: Data and code for the article "Experimental addition of marine-derived nutrients affects wildflower traits in a coastal meta-ecosystem" by Allison M. Dennert, Elizabeth Elle, and John D. Reynolds
---

This dataset is associated with a long-term field experiment conducted in Haíɫzaqv Territory on the Central Coast of British Columbia, Canada. The experiment has a randomized block design with four treatments, and 25 replicates per treatment. Wildflower leaf isotopes, leaf area, floral display size, and seed production were quantified over three summer growing seasons (2017-2019) in order to understand the relationship between marine-derived nutrients from salmon and drift seaweed and plant growth and reproduction. Up to 3 individuals of each focal plant species were measured per treatment plot. Individuals are nested within plot, nested within block, except for the isotope data where one individual of each species was sampled


## Description of the Data and file structure

There are two data files included: 2017leafisotopes.csv, and focaltraits.csv

Details for: 2017leafisotopes.csv
--------------------

* Description: A comma-delimited file containing the results of stable isotope analyses for carbon and nitrogen in plant leaves collected from each experimental treatment plot. One leaf of each of four plant species was collected from a subset of the treatment plots in each of the three years that the experiment was conducted.

* Format(s): .csv

* Missing Data Codes: These data contain no missing values. 

* Dimensions: 630 rows x 10 columns

* Variables: 
	* year: The year (YYYY) that each leaf sample was collected
	* block: The experimental block (possible values 1-25) that each treatment plot belongs to
	* treatment: The experimental treatment applied to a 1 x 1 m plot; A = pink salmon (Oncorhynchus gorbuscha) carcass, B = drift seaweed (Fucus distichus), C = both salmon carcass and seaweed, D = control
	* species: A four-letter species code denoting the plant species that the leaf was collected from; ASSU = Aster subspicatus (synonym: Symphyotrichum subspicatum), ACMI = Achillea millefolium, CAMI = Castilleja miniata, POAN = Potentilla anserina
	* sample.id: An amalgamated variable comprised of the specimen's collection year, block, treatment, and species in a 13-character ID (YYYY-species-blockplot)
	* d13C: The amount of the carbon-13 isotope present in the leaf sample, presented in parts per thousand or units per mil (‰)
	* C.ug: The total mount of carbon present in the leaf sample, presented in micrograms (μg)
	* d15N: The amount of the nitrogen-15 isotope present in the leaf sample, presented in parts per thousand or units per mil (‰)
	* N.ug: The total mount of nitrogen present in the leaf sample, presented in micrograms (μg)
	* sample.mg: The weight of the sample in milligrams (mg)
	

Details for: focaltraits.csv
--------------------

* Description: A comma-delimited file containing the measured values of plant and floral traits in wildflowers growing in the experimental treatment plots. Up to three individuals of each focal plant species was tagged in each treatment plot and measured throughout the growing season in each of the three years that the experiment was conducted. 

* Format(s): .csv

* Missing Data Codes: missing data is coded as "NA"

* Dimensions: 2196 rows x 34 columns

* Variables: 
	* id: An amalgamated variable comprised of the specimen's sampling year, block, treatment, and species in a 13-character ID (YYYY-blockplot-species-individual)
	* year: The year (YYYY) that each leaf sample was collected
	* date: The date (dd-Month-yy) that an individual's leaf and floral display size was measured
	* block: The experimental block (possible values 1-25) that each treatment plot belongs to
	* plot: The experimental treatment applied to a 1 x 1 m plot; A = pink salmon (Oncorhynchus gorbuscha) carcass, B = drift seaweed (Fucus distichus), C = both salmon carcass and seaweed, D = control
	* plot.id: An amalgamated variable comprised of the specimen's block and plot (blockplot)
	* species: A four-letter species code denoting the plant species that the leaf was collected from; ASSU = Aster subspicatus (synonym: Symphyotrichum subspicatum), ACMI = Achillea millefolium, CAMI = Castilleja miniata, POAN = Potentilla anserina
	* indiv: The individual number assigned to each plant within a plot (up to 3 individuals of each species were tagged in each plot)
	* leaf. length: the length of a randomly selected leaf (mm)
	* leaf.width: the widest width of a randomly selected leaf (mm)
	* total.leaves: the total number of leaves on a plant
	* flower.number: the total number of flowers on a plant
	* display.size.1-8: the floral display size of each inflorescence on a plant (mm); for individuals with more than one inflorescence (8 maximum), each successive flower was measured such that individuals in the data set with missing values did not grow that many flowers
	* seed.set.1-14: the number of seeds produced by each inflorescence on a plant (mm); for individuals with more than one inflorescence (14 maximum), each successive flower's seeds were counted such that individuals in the data set with missing values did not grow that many flowers
	

## Sharing/access Information

Links to other publicly accessible locations of the data: https://github.com/adennert/wildflower-traits

These data were generated by the study's co-authors and not derived from another source.







