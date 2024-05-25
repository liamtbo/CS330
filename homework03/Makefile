cc = gcc

sources = main.c mmio.c
# what the -o flag compiles the c files into.
objects = $(addsuffix .o, $(basename $(sources)))

flags = -g -std=c11 

# executable name -> write ./spmv <args> to run program
target = spmv

all: $(target) 
# target = spmv which  is the executable file
# compiles object files, .o, into executable
# objects -> addsuffix .o -> %.o creates object file
$(target) : $(objects)
	$(cc) $(flags) -o $(target) $(objects)

# whenever a .o is seen with anything, %, infront of it, this gets
# invoked
%.o : %.c
	$(cc) -c $(flags) $< -o $@

clean:
	rm -rf $(target) $(objects) 
