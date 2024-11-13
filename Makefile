# Compiler and flags
MPICXX = mpicxx
CXXFLAGS = -O2 -std=c++11 -Wall -fopenmp

# Executable name
TARGET = spellcheck

# Source files
SRCS = a2_v3.cpp hash_table_lock.cpp  
# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJS)
	$(MPICXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Compile source files into object files
%.o: %.cpp  
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)



.PHONY: all clean run
