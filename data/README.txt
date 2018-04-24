- MAX_Cell#_signal_reduction.csv:
These are ImageJ-generated files containing information on how the green channel decreases through time for each large image.

- mitophagy_summary_intensity_mean_ch2.csv:  ##### IMPORTANT ######
This data set combines ch2_mean_intensity from the results in mitophagy_results_from_imaris/.

- Synchronisation
These mitophagy time courses tend to have an peakillatory behaviour but also include a fair amount of noise. 
(1) In order to attempt a synchronisation, these time courses were splined smoothly (spar=0.4) so that the main time course profile was more identifiable (FILE: mitophagy_summary_intensity_mean_ch2_spline.csv). 
(2) After computing the splines, the synchronisation was done manually, highlighting the time regions of the peaks (FILE: mitophagy_summary_intensity_mean_ch2_spline__manual_synchronisation.ods). 
(3) The table of synchronising the splined time courses was then replaced with the original table of time courses (non splined!) (FILE: mitophagy_summary_intensity_mean_ch2__manual_synchronisation.ods). 
(4) The tables of splined and original synchronised time courses were then used for plotting and further analysis (FILES: mitophagy_summary_intensity_mean_ch2__synchronised.csv, mitophagy_summary_intensity_mean_ch2_spline__synchronised.csv).
(5) The FILES: mitophagy_summary_intensity_mean_ch2__synchronised_filtered.csv, mitophagy_summary_intensity_mean_ch2_spline__synchronised_filtered.csv contain the following edits (via R script: synchronised_time_courses.R) 
    (a) the time series 3A, 3B, 5A 6D, 6F2, and 16sA were removed because too irregular;
    (b) time course cut off at [0, 520]s and [1640,1720]s because there was no information in there;    
    (c) the time series were log10-transformed;
    (d) min-max rescaling was applied for each time course;

- the files with "_regularised" contain the data where the green channel does not decrease over time. These data were regularised using the mean regression line in folder 6_synchronised_time_courses.

- selected_movies__17.csv contains the name of the selected frames and the mitochondrial diameter data.
