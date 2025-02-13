#!/usr/bin python3

# Import the necessary libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
from matplotlib import gridspec

# library &amp; dataset
import seaborn as sns


def main():

	sns.set_context("paper", font_scale=1.5)

	filename = sys.argv[1]
	dataset_name = sys.argv[2]
	plen = sys.argv[3]
	infile_name = filename
	#outfile_name = filename + "_Time_ecoli.png"
	outfile_name = filename + "_time_locate_" + dataset_name + ".png"

	# Load in data
	df = pd.read_csv(infile_name)

	# Initialize Figure and Axes object
	fig = plt.figure(figsize=(6, 5))
	ax = fig.add_subplot()

	#The cubehelix color palette system makes sequential palettes with a linear increase or decrease in brightness and some variation in hue.
	#This means that the information in your colormap will be preserved when converted to black and white (for printing) or when viewed by a colorblind individual.
	# sns.set_palette("colorblind")
	sns.set_palette("bright")

	df = df[ df["Dataset"] == dataset_name ]
	df = df[ df["pattern length"] == int(plen) ]

	#g = sns.lineplot(x="space", y="time", hue="index", style="pattern length", markers=True, dashes=True, data=df, ax=ax, **{'markersize':10}).set_title("P&C cere")
	g = sns.lineplot(x="space", y="time", hue="index", style="index", markers=True, dashes=True, data=df, ax=ax, **{'markersize':10}).set_title("P&C cere - pattern length "+str(plen))

	ax.set_yscale('log', base=10)
	#ax.set_xscale('log', base=2)
	ax.grid(which='major', axis='both')
	ax.set_ylabel("Time to locate 100k patterns [sec]")
	ax.set_xlabel("Resident set size [GB]")
	ax.set_xlim(-100000000,2300000000)
	plt.minorticks_off()

	# Removing "Datasets" from legend
	handles, labels = ax.get_legend_handles_labels()
	#ax.legend(title="Sigma", title_fontsize='11', handles=handles[0:], labels=labels[0:], loc='upper right', fontsize='10')
	'''
	labels[1] = "$r$-index"
	labels[2] = "extended $r$-index"
	labels[3] = "$sr$-index"
	labels[4] = "CHICO"
	labels[5] = "grammar index"
	'''
	ax.legend(handles=handles[0:], labels=labels[0:], loc='upper right', fontsize='9')

	fig.tight_layout()  # otherwise the right y-label is slightly clipped

 	# Save the figure
	plt.savefig(outfile_name, dpi=300)#,  additional_artists=art, bbox_inches="tight")#bbox_extra_artists=(lgd,), bbox_inches='tight')
	#xplt.savefig(outfile_name + ".pdf", bbox_inches='tight')


if __name__ == '__main__':
	main()
