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
	Parameters.java \
	IOUtil.java \
	Constants.java \
	ForceModel.java \
	SPMModel.java \
	NBIModel.java \
	MONDModel.java \
	Integrator.java \
	MathUtil.java \
	SetupUtil.java \
	QuaternionRotator.java \
	RotationController.java \
	ScatterPanel.java \
	Run.java

default: all

all: classes jar

initdirs:
	mkdir -p build
	mkdir -p build/classes
	mkdir -p build/jar

classes: initdirs $(CLASSES:.java=.class)

clean:
	$(RM) build

jar:
	jar cvf build/jar/jspamga.jar -C $(CDIR) edu 
