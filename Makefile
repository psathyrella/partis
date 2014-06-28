VERSION ?= 0.1.0
ARCH ?= amd64

PYSAM_VERSION ?= 0.7.8
NETWORKX_VERSION ?= 1.9

GIT_DESCRIBE = $(shell git describe --always)

package: pysam networkx ighutil vdjalign

upload: package
	(cd debs && \
		bundle exec deb-s3 upload --sign 3AA09EC1 -a amd64  *.deb --bucket cmccoy-debian-repo)

pysam: debs/python-pysam_$(PYSAM_VERSION)_$(ARCH).deb

networkx: debs/python-networkx_$(NETWORKX_VERSION)_all.deb

ighutil: debs/ighutil_$(VERSION)_$(ARCH).deb

vdjalign: debs/vdjalign_$(VERSION)_$(ARCH).deb

debs/ighutil_$(VERSION)_$(ARCH).deb:
	mkdir -p debs
	+make -C clj all
	mkdir -p tmp/usr/local/bin
	cp clj/bin/ighutil tmp/usr/local/bin
	@echo "Packaging $@"
	(cd debs && bundle exec fpm \
		-C ../tmp \
		--version $(VERSION) \
		-d 'java7-runtime-headless' \
		-n 'ighutil' \
		--license GPLv3 \
		--maintainer "Connor McCoy <cmccoy@fhcrc.org>" \
		--description "Tools for working with output of vdjalign. sha: $(GIT_DESCRIBE)" \
		--url "http://github.com/cmccoy/ighutil" \
		--deb-suggests vdjalign \
		--category universe/math \
		-s dir \
		-t deb \
		.)
	rm -rf tmp

debs/python-networkx_$(NETWORKX_VERSION)_all.deb:
	mkdir -p debs
	@echo "Packaging $@"
	(cd debs && bundle exec fpm \
		-s python \
		-t deb \
		networkx)

debs/python-pysam_$(PYSAM_VERSION)_$(ARCH).deb:
	mkdir -p debs
	@echo "Packaging $@"
	(cd debs && bundle exec fpm \
		--no-python-dependencies \
		-d 'cython (>=0.19)' \
		-s python \
		-t deb \
		pysam)

debs/vdjalign_$(VERSION)_$(ARCH).deb:
	mkdir -p debs
	@echo "Packaging $@"
	(cd debs && bundle exec fpm \
		-n vdjalign \
		-s python \
		-t deb \
		-d 'samtools' \
		--license GPLv3 \
		--maintainer "Connor McCoy <cmccoy@fhcrc.org>" \
		--description "Align IGHV gene sequences. sha: $(GIT_DESCRIBE)" \
		--url "http://github.com/cmccoy/ighutil" \
		--category universe/math \
		--deb-suggests ighutil \
		../python/setup.py)

docker-test: Dockerfile package
	docker build .

.PHONY: package docker-test
