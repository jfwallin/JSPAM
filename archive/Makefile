# example Makefile from 
# http://myweb.stedwards.src/src/edu/laurab/help/javamakefile.html
#JFLAGS = -g
JFLAGS = 
JC = javac
CDIR = build/classes
SDIR = ./ 
CP = .:build/classes
RM = /bin/rm -rf
.SUFFIXES: .java .class
.java.class:
	$(JC) $(JFLAGS) -sourcepath $(SDIR) -classpath $(CP) -d $(CDIR) $*.java

CLASSES = \
	src/OrbitUtil.java \
	src/JSpamParams.java \
	src/JSpam.java \
	src/DFJSpam.java \
	src/Driver.java

default: all

all: classes jar

initdirs:
	mkdir -p build
	mkdir -p build/classes
	mkdir -p build/jar

classes: initdirs $(CLASSES:.java=.class)

clean:
	$(RM) build
	$(RM) output 

jar:
	jar cvf build/jar/jspamarchive.jar -C $(CDIR) edu 
