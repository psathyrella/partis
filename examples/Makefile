all: pf-example rare-events-example cluster-example

.PHONY: pf-example rare-events-example cluster-example clean

clean:
	$(MAKE) -Cpf clean
	$(MAKE) -Crare-events clean
	$(MAKE) -Ccluster clean

pf-example:
	$(MAKE) -Cpf all

rare-events-example:
	$(MAKE) -Crare-events all

cluster-example:
	$(MAKE) -Ccluster all
