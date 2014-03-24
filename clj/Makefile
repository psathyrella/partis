# ighutil -- Create executable script

all:
	rm -f target/*.jar
	lein uberjar
	cat bin/ighutil.template target/*-standalone.jar > bin/ighutil
	chmod +x bin/ighutil

.PHONY: all
