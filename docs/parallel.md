[Up to table of contents](contents.md)

### Parallelization

###### In general

The number of processes on your local machine is set with `--n-procs N`, which defaults to the number of cpus.

In order to parallelize over more processes than the local machine can handle, we currently support slurm and sge: specify one or the other with `--batch-system`.
The default options for each should work, but if you need to add extras (for instance to reserve particular memory requirements) use `--batch-options`, e.g. `--batch-options="--foo bar"`.
Note that you should not specify stdout/stderr locations (`-e` or `-o`) for sge -- we add these options automatically because we need to be able to parse them for each job.
By default, partis writes its temporary files to a working directory of the form `/tmp/$USER/hmms/$RANDOM`.
If you're running on a batch system, though, you need the working directory to be a network mount that every node can see; set this with `--workdir`, e.g. `--workdir /path/to/nfs/$USER/hmms/$RANDOM`.

###### annotate

Sequence annotation lends itself quite readily to independent parallelization.
It should take 0.1-1 second per sequence; if you want it to go faster, just increase `--n-procs`.

###### partition

Clustering, however, really doesn't nearly as easily lend itself to independent parallelization -- we need to, at least approximately, compare each sequence to every other one.
For the full partis method, we get around this by starting with `--n-procs` processes.
The input sequences are split evenly among these, and each process does all-against-all comparison (with many optimizations to avoid the full likelihood calculation) of all of its allotted sequences.
The results of this first round are collected and merged together, and then reapportioned among a new, smaller, number of processes.
This is continued until we arrive at one final process which is comparing all sequences.
Since at each stage we cache every calculated log probability, while the later steps have more sequences to compare, they also have more cached numbers at their disposal, and so it's possible to make each step take about the same amount of time.
We currently reduce the number of processes by about 1.6 at each step, as long as the previous step didn't have to calculate too many numbers.

For the vsearch partis method (`--fast`), vsearch does all its usual cleverness to avoid all-against-all comparison, and is thus blindingly fast.
