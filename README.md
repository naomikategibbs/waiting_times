# waiting_times
This is the R code for a decision analytic model to estimate the health impact of waiting for elective surgery in England.
This version of the model is deterministic. It draws in a number of parameter inputs estimated using routine datasets such 
as Hospital Episode Statistics (HES), Office for National Statistics mortality figures, Patient Reported Outcome Measures (PROMs) and 
Private Health Information Network. It also draws on data from the academic literature.

Only the model is provided, and the data is not shared as some of the model inputs are estimated from confidential, protected data, namely the HES and PROMs data.

The model is primarily a Markov model. 

The full methodological approach and all the data inputs are outlined in the following paper: "include link to published paper."

For further information please get in touch with naomi.gibbs@york.ac.uk

This code was developed by Naomi Gibbs with support from Simon Walker and Susan Griffin. 
The underlying Markov model utilised code from this repository https://github.com/LSHTM-GHECO/DM4HEE_RCode, originally created to write R code that was complementary to the "Decision Modelling for Health Economic Evaluation" course. See this website for further information about the course: https://modelling.he-evalcourses.com/
It was briefly reviewed by Naomi Gibbs and Julia Hatamyar. It has been tested for extreme values and across eight different elective procedures.

This work is published under a CC-BY license
