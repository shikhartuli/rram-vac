import sys
import os 


ecg_input = open("ecg_signal","r")


line = ' '

output_counter = 0
files_counter = 0

while files_counter < 20:
	outputname = "ecg_3s_" + str(files_counter)
	ecg_output = open(outputname, "w")
	while output_counter < 1500:
		line = ecg_input.readline()
		ecg_output.write(line)
		output_counter = output_counter + 1
	ecg_output.close()
	cmd =  './compress ' + outputname + ' out.csv >> report_global_' + str(files_counter) +'.csv'
	os.system(cmd)

	files_counter = files_counter + 1
	output_counter = 0
