This folder contains real mitophagy events which were analysed by Piero. This data set is used for further analysis.




- data/:
data sets for data analysis.

- 1_large_images/:
confocal images AFTER colour alignment using imagej. 

- 2_selected/:
mitophagy events selected from movies in large_images/

- 3_mitophagy_results_from_imaris/: 
Imaris generated statistics after tracking and quantification.

- 4_ signal_decrease_data/:
time course of green signal decrease. These represent how the green channel signal decreases over time. Generated with imagej.

- 5_mitophagy_time_courses_plots/:   ##### IMPORTANT FOLDER #####
Plot each mitophagy event time course (ch2_mean_intensity, data: data/mitophagy_summary_intensity_mean_ch2.csv). 
Calculate and plot the spline (spar=0.4) for each mitophagy event time course (ch2_mean_intensity, generated data: data/mitophagy_summary_intensity_mean_ch2_spline.csv)

- 6_synchronise_time_courses/:
Synchronise and plot time course data. 

- 7_time_courses_w_adjust_green_ch/:
Adjust the time courses so that there is no signal decline over time. This was not background signal. 


- movies_orthogonal_atg13_circle.csv:
File containing the list of selected movies showing an accurate circle.




Others:

- movies_with_tracking_on/:
some movies generated with Imaris after tracking.

- test_imagej_quantif/: 
a quick green channel signal quantification using ImageJ. This is standalone and independent from the other analyses.
