## ⚠️ **CRITICAL: Virtual Environment Required**

**ALWAYS activate the virtual environment before running any commands:**
```bash
source .venv/bin/activate
```

## docs

Read the `docs/` dir to understand how the code base works, and how to install it.
Usually, install from source with `pip -e` for editable install.


## testing

See docs for testing details.
The minimal test is:

partis-test.py --quick

Any significant changes require the standard test that runs more actions:

    partis-test.py

And any changes that affect paired code should instead run (generally don't need both non-paired and paired):

    partis-test.py --paired

Any crashes obviously need to be fixed.
The color-coded output tells you if either results or time required have changed: red is a big change, yellow is a smaller change.
Right now there's a few differences that need to be updated in the test framework (e.g. time required is larger than in the ref results because the ref times required need to be updated for changing a default option).
Some files also differ in the parameter and simulation dirs, and need updating.


## notes

Do **NOT** remove comments, commented code, or TODOs unless absolutely certain that you are making changes that fix them or make them irrelevant.
Comments and TODOs are not 'cruft' to remove, they are purposefully placed to remind of things in the future.

Also do not remove functionality or checks that seem peripheral without CAREFULLY asking whether they should be removed.

## data

The parameter and simulation dirs under `test/` (`test/ref-results/`, `test/new-results/`, and their
`-slow` counterparts) are **test fixtures, not data**. They are tiny -- e.g.
`test/ref-results/test/parameters/data/hmm` is built from 47 sequences -- and they are not
representative of anything. Never draw a conclusion about partis's behaviour on real repertoires from
them, and never quote a number measured on them as if it described data. They exist to detect changes
in output, and that is all they are good for. Measure on real repertoires instead.

## judging changes

BCR repertoires vary along many axes at once: mutation rate, clonal family sizes, tree shape, germline
set, allele frequencies, read length, locus, sample size. That space can't be scanned, even in
principle -- there is no way to enumerate tree shapes, for instance. So "it worked on a few samples"
does not by itself justify a change. Before making one, work out from how the code works which inputs it
affects and how, including the ones that aren't at hand (very high and very low mutation, tiny and
huge samples, rare alleles, other loci), and make sure it doesn't break some other part of that space.
Samples can confirm that understanding. They can't replace it.
