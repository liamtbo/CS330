# all the .c files in the current directory
sources = $(wildcard *.c)
# adds an .o to all of these files
objects = $(addsuffix .o, $(basename $(sources)))
flags = -g -W -Wall -std=c11 -DDEBUG=0 -DTEST=0
# bfs is name of final executable
target = bfs 


$(target) : $(objects)
	gcc -o $(target) $(objects)

# if make encounters a target dependency "example.o", and there exists
# a corresponding source file example.c, this will be called
%.o : %.c
	gcc -c $(flags) $< -o $@
clean :
	rm $(target) $(objects)

run : bfs test2/test2.txt
	./bfs test2/test2.txt 75