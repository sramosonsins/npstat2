all: npstat

npstat: 
	gcc -o npstat2 NPStat-v2.c -lgsl -lgslcblas -lm -lhts
clean: 
	rm npstat
