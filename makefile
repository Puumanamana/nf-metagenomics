CONF ?= local,docker,test

test:
	nextflow main.nf -resume -profile $(CONF)
clean:
	rm -rf work .nextflow* output_test
