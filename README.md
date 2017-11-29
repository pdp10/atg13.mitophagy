
# ATG13 mitophagy data analysis


### Summary
Data analysis for mitophagy events targetting ATG13 as readout. These data were quantified and analysed by me (Piero).


### Pipeline structure

- data/: data sets for data analysis.

- 1_large_images/: confocal images AFTER colour alignment using imagej. 

- 2_selected/: mitophagy events selected from movies in large_images/

- 3_mitophagy_results_from_imaris/: Imaris generated statistics after tracking and quantification.

- 4_ signal_decrease_data/: time course of green signal decrease. These represent how the green channel signal decreases over time. Generated with imagej.

- 5_mitophagy_time_courses_plots/: plot each mitophagy event time course (data: data/mitophagy_summary_intensity_mean_ch2.csv). Calculate and plot the spline (spar=0.4) for each mitophagy event time course (generated data: data/mitophagy_summary_intensity_mean_ch2_spline.csv)

- 6_synchronise_time_courses/: synchronise and plot time course data. 

- 7_time_courses_w_adjust_green_ch/: adjust time courses so that there is no signal decline over time. Note: this is not background signal. 

- 8_delay_analysis/: extract delay times considering the highest (peaks) and lowest values of the mean of the oscillations.

- 8_upper_lower_bound_analysis/: extract the upper and lower values of the oscillations.

- 9_time_courses_data_for_copasi/: generate the data sets for parameter estimation using Copasi.

- 10_tc_heatmap/: plot time courses in a heatmap.


### Other folders

- movies_with_tracking_on/: some movies generated with Imaris after tracking.

- test_imagej_quantif/: a quick green channel signal quantification using ImageJ. This is standalone and independent from the other analyses.

