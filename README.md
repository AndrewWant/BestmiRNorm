# BestmiRNorm
## Novel method for identifying miRNA normalisers for qPCR
This script enables more transparent identification of the optimal normalisers for a qPCR experiment.
It accepts an input file, with a pre-specified format, for the qPCR data along with information about the potential normalisers in the dataset. 
The current implementation can examine up to 11 normalisers in a time-efficient manner. If more than 11 normalisers are analysed, this will take at least 45 minutes, and more than 15 will consume considerable system resources.

This repository also contains an example Microsoft Excel (.xlsx) file for reference.

This code may be adapted in the future to reduce this issue by batching the processing.

The current implementation considers the following:
*  Kolmogorov-Smirnov score, evaluating distribution similarity of two experimental/disease classifications
*  Average absolute distance of class means from 0 
*  Average of class standard deviations

And these scores are calculated for all possible combinations of the potential normalisers. All normaliser combinations are then ranked, with the lowest value assigned to the combination with the highest value for a particular scoring criterion. These three ranked values may then be multiplied by a weighting factor which represents the experimenter's assessment of the importance of that criterion to the overall assessment of stability. This system was validated with integer values between 1 and 3, although in theory any combination will work. 

Plotting tools may be made available at a future time, depending on specific requests. Crucially, additional functionality can be implemented by anyone due to the implementation in a widely used open-source programming language. 
