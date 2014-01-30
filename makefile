#makefile in cuvodec
NAME= decoders
BIN_DIR= ./bin
SRC_DIR= ./src
EXE_DIR= ./exe
INCLUDE_DIR= ./include
MAKE_DIR= ./make
INCLUDE= -I$(INCLUDE_DIR)

TARGET= $(EXE_DIR)/TestMassDecoder

vpath %.cpp $(SRC_DIR)
vpath %.cu $(SRC_DIR) 
vpath %.h $(INCLUDE_DIR)
vpath %.d $(MAKE_DIR)


CPP_OBJS= $(patsubst $(SRC_DIR)/%.cpp, $(BIN_DIR)/%.o, $(wildcard $(SRC_DIR)/*.cpp))
CU_OBJS= $(patsubst $(SRC_DIR)/%.cu, $(BIN_DIR)/%.o, $(wildcard $(SRC_DIR)/*.cu))

CPP_MAKE_TARGETS= $(patsubst $(SRC_DIR)/%.cpp, $(MAKE_DIR)/%.d, $(wildcard $(SRC_DIR)/*.cpp))
CU_MAKE_TARGETS= $(patsubst $(SRC_DIR)/%.cu, $(MAKE_DIR)/%.d, $(wildcard $(SRC_DIR)/*.cpp))

TARGET_LIB= $(BIN_DIR)/$(NAME)_lib.a

.PHONY: link compile clean document makegen

link: compile $(TARGET)
	@echo -- Link Ended

clean:
	rm -f $(BIN_DIR)/*.o
	rm -f $(MAKE_DIR)/*
	rm -f $(EXE_DIR)/*
	@echo "--> Clean Ended"

document:
	@echo "--> Documentation to be done" 

compile: $(CPP_OBJS) $(CU_OBJS) $(CPP_MAKE_TARGETS) $(CU_MAKE_TARGETS)
	@echo "Compiling ended for.... " $(CPP_OBJS) $(CU_OBJS)

#include $(CPP_MAKE_TARGETS)
#include $(CU_MAKE_TARGETS)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "Compiling... $@"
	g++ -Wall -c $< $(INCLUDE) -o $@ 
	
$(BIN_DIR)/%.o: $(SRC_DIR)/%.cu
	@echo "Compiling... $@"
	nvcc -c $< $(INCLUDE) -o $@
	
$(TARGET): $(TARGET_LIB) $(CU_OBJS)
	g++ -Wall $^ $(INCLUDE) -o $@ -lcudart
	
$(TARGET_LIB): $(CPP_OBJS)
	ar rvs  $(TARGET_LIB) $(CPP_OBJS)

makegen: $(CPP_MAKE_TARGETS) $(CU_MAKE_TARGETS)
	@echo "Making ended for ..." $(CPP_MAKE_TARGETS) $(CU_MAKE_TARGETS)

$(MAKE_DIR)/%.d: $(SRC_DIR)/%.cpp
	@echo "Making..... $@" 
	@set -e;\
	rm -f $@;\
	g++ -MM $(INCLUDE) $< >$@.$$$$;\
	sed 's,\($*\)\.o[ :]*, $(BIN)/\1.o $@: ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$

$(MAKE_DIR)/%.d: $(SRC_DIR)/%.cu
	@echo "Making..... $@" 
	@set -e;\
	rm -f $@;\
	nvcc -MM $(INCLUDE) $< >$@.$$$$;\
	sed 's,\($*\)\.o[ :]*, $(BIN)/\1.o $@: ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$
