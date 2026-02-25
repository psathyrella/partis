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

**NEVER* stage or commit files to git. Staging and committing is where I check your changes. You can use git diff and status commands as much as you like.
