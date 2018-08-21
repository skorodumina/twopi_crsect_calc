
CXX           = g++ -Wno-write-strings -Wno-pragmas
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINCLUDE  := -I$(shell root-config --incdir)



all:
	make script_exe

script_exe: script_xsect_calc.$(ObjSuf) 
	$(CXX) -g -o script_exe $^ $(ROOTGLIBS)  



%.$(ObjSuf): %.$(SrcSuf)
	$(CXX) -g -c $(ROOTINCLUDE) -c $<

clean:
	rm -f *.o
	rm -f *Dict.*
	rm -f G__*
	rm script_exe
	

