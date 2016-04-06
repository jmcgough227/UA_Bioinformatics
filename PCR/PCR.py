import os
import subprocess

sequence_file = open("sequence.txt", "r")
rna = sequence_file.read()
sequence_file.close()
rna = rna.upper()

start_pos = raw_input('Starting position of the desired strand --> ')

while 1:
	if start_pos.isdigit():
		if int(start_pos) <= len(rna) - 200:
			break;
	start_pos = raw_input('Invalid entry. Please enter a positive integer between 0 and %d --> ' % (len(rna) - 200))

start_pos = int(start_pos)


primer_config_file = open("primer_config.txt", "w")

primer_config_file.write("SEQUENCE_ID=ebola\n")
primer_config_file.write("SEQUENCE_TEMPLATE=" + rna + "\n")
primer_config_file.write("SEQUENCE_INCLUDED_REGION=%d,%d\n" %(start_pos, start_pos + 200))
primer_config_file.write("PRIMER_NUM_RETURN=1\n")
primer_config_file.write("PRIMER_TASK=generic\n")
primer_config_file.write("PRIMER_PICK_LEFT_PRIMER=1\n")
primer_config_file.write("PRIMER_PICK_RIGHT_PRIMER=1\n")
primer_config_file.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + os.getcwd() + "/primer3-2.3.6/src/primer3_config/\n=")

primer_config_file.close()

subprocess.call(['primer3-2.3.6/src/primer3_core', '-output='+os.getcwd()+'/primer3_output.txt', 'primer_config.txt'])

primer_output_file = open('primer3_output.txt', 'r')

left_primer = ""
right_primer = ""

for line in primer_output_file:
	if 'PRIMER_LEFT_0_SEQUENCE=' in line:
		left_primer = line.replace('PRIMER_LEFT_0_SEQUENCE=','')
		left_primer = left_primer.replace('\n','')
	elif 'PRIMER_RIGHT_0_SEQUENCE=' in line:
		right_primer = line.replace('PRIMER_RIGHT_0_SEQUENCE=','')
		right_primer = right_primer.replace('\n','')

primer_output_file.close()

print "left primer: " + left_primer
print "right primer: " + right_primer

