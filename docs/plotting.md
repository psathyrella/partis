#### plotting

Note that in order to make plots for the `partition` action, you have to have R installed, along with several [extra packages](install.md#plotting).

The addition of `--plotdir <plotdir>` to most partis commands will write to disk a variety of plots that we find useful for that command.
These plots are written as svg files to various subdirectories of `<plotdir>`, along with html files displaying clickable summaries of the svgs.
You typically want to view the html files in a browser, so a good way to see what's available might be to run `find <plotdir> -name '*.html' | xargs firefox`.
In addition, plots that are simply histograms usually also have their histogram content written to a csv in the same directory.
This makes it easier for later comparisons across several directories (for example with `bin/compare-plotdirs.py`).

Adding `--plotdir <plotdir>` to `cache-parameters` or `annotate` will write plots summarizing the distributions of rearrangement parameters in the sample (gene calls, deletion and insertion lengths, SHM properties, etc.).
Examples of these plots can be found [here](http://psathyrella.github.io/partis/example-plots/hmm/overall.html) for repertoire-wide rearrangement parameters, and [here](http://psathyrella.github.io/partis/example-plots/hmm/mute-freqs/overall.html) for repertoire-wide SHM frequencies broken down in various ways.
In order to view similar files that you've written to your file system, you'd typically run, for example, `firefox docs/example-plots/*/*.html` (or `google-chrome`), where we've included copies of the same plots into a subdirectory of this documentation.
If you set `--make-per-gene-plots` in addition to `--plotdir`, similar plots will also be written breaking things down for each individual gene, for instance the per-position SHM rates and 5' deletion lengths for each V gene.

An overview of the plots written for the `partition` action can be found [here](http://psathyrella.github.io/partis/example-plots/partitions/overview.html).
At top left is a plot showing a colored blob/slug for each clonal family (sorted by size) whose extent along the x direction shows the distribution of SHM rates within that family.
The (zero-indexed) family rank and size are shown along the right side.
Note that the three colors (green, blue, and yellow) have no separate significance, and are used only to distinguish adjacent slugs.
This particular plot also shows the result of setting some sequences of interest using `--queries-to-include a:b:z:ac`, such that sequences labeled a, b, z, and ac will be highlighted in red.
A typical use case for this option is if you have made several previous seed partition runs (with e.g. `--seed-unique-id a`), and you want to see how the families of the seed sequences fit into the larger repertoire.
Only the first few of these slug plots (with the biggest clusters) are shown in `overview.html` -- the rest are in the `shm-vs-size/` subdirectory (or click `shm-vs-size` link at top of page).
The rest of the top row is occupied by log and non-log scatter plots showing the size vs mean SHM rate for each family.

Below this, there is a "multi-dimensional scaling" (MDS) plot for each clonal family, where each sequence in each family is a point on that family's plot.
For a good description of MDS, you should go to the google, but it's essentially an improved version of principal component analysis (PCA).
For our purposes, MDS takes each family (as a cluster in 400-odd dimensional sequence/hamming distance space) and squishes it out into two dimensions, choosing axes such as to maximize how squished out the family gets, while as nearly as possible preserving each inter-sequence distance from the real, 400-odd dimensional space.
This means that while there is no easy biological interpretation for the various directions on these plots, they do a good job of giving, at a glance, an overview of the basic structure of each family.
The inferred naive sequence for each cluster is shown as a red point on these plots (as are any queries specified with `--queries-to-include`), and the SHM rate of each sequence is given by its transparency (darker is more mutated).
So in cases where dots get uniformly less transparent as they get further from the red naive dot, this tells you that the dimension reduction is not losing very much information.
In cases where the plots are, on the other hand, uniformly speckled all over, the sequences are distributed more evenly across the 400-odd dimensional space (i.e. there wasn't a way to squish down to two dimensions without losing lots of information).
The plot title shows the family's size (this size, and the plot's file name, match the numbers on the right side of the slug plots), as well as the amino acid translation of its CDR3.
The overview html again only shows plots for the largest few clusters, while the rest can be found in the `mds/` subdirectory (or click `mds` link at top of page).

For a description of the plots written during germline inference, see [here](germline-inference.md).
