localrules: all, testpipe

TARGETS = ('/ei/scratch/linsmith/TESTASSEM/out.touch')

rule all:
	input: TARGETS

rule testpipe:
        input:
		"/ei/scratch/linsmith/TESTASSEM/in.touch"
        output:
                "/ei/scratch/linsmith/TESTASSEM/out.touch"
        threads:
                1
        shell:
                "touch {output}"
