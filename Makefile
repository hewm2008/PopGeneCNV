CXX=g++
BIN := ./bin
SRC := ./src

LDFLAGS=-lz  -lpthread -llzma	 -lbz2  -lcurl	 -lncurses	 -lhts 
INCLUDE=-I./include -I/hwfssz4/BC_PUB/Software/08.Centos7/htslib-1.16/include/
LIBS=-L./lib  -L/hwfssz4/BC_PUB/Software/08.Centos7/htslib-1.16/lib/

CXXFLAGS= -g  -O3  $(INCLUDE) $(LIBS)
#  you can add the [ -I /path/.../samtools-1.6/htslib-1.6   -L /path/.../samtools-1.6/htslib-1.6 ] here #no need by fqcheck


all: $(BIN)/PopGeneCNV

$(BIN)/PopGeneCNV: $(SRC)/PopGeneCNV.o 
	$(CXX)   $^ -o $@   $(LDFLAGS)  $(INCLUDE) $(LIBS)

$(SRC)/%.o: %.cpp
	$(CXX)  -c $(CXXFLAGS)  $< -o $@ 

clean:
	$(RM) $(BIN)/*.o  $(SRC)/*.o

