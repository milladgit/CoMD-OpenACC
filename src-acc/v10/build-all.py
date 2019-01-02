
import os, sys

MIN=512
MAX=512

for i in range(MIN, MAX+1):
	print "Building for maxregcount %d ... " % (i)
	os.system("make -f Makefile.tesla clean &> /dev/null")
	os.system("MAXREGCOUNT=%d make -f Makefile.tesla &> build.output" % (i))
	os.system("grep -in \"ptxas warning : Too big maxrregcount value specified\" build.output > building-%d.result" % (i))
	print "Building for maxregcount %d ... Done --- See output file: building-%d.result" % (i, i)


