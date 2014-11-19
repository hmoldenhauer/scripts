all: file_list.txt ampl.txt

file_list.txt:
	ls *.dat >> file_list.txt

ampl.txt: modelocking_fit_ampl_listing.py
	python modelocking_fit_ampl_listing.py
