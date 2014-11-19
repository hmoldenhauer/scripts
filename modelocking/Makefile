all: file_list.txt ampl.txt

file_list.txt:
	ls data/*.dat >> file_list.txt

ampl.txt: modelocking_fit_ampl_listing.py
	mkdir -p fit
	python modelocking_fit_ampl_listing.py
