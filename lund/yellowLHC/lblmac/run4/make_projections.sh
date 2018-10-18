#!/bin/bash

for radius in 0.2 0.3 0.4 0.5
do
	./make_lund_plots.py merged.root -r ${radius}
	./make_lund_plots.py merged.root --charged -r ${radius}
done
