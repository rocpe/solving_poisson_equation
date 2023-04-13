# *****************************************************
# compiler
CC = g++

CFLAGS  = -g -Wall -pedantic -larmadillo -O3 

# The build target 
TARGET = a.out
RUN = ./
R = Rscript
RSCRIPT = plots.r

# ****************************************************
 
main: main.o 
	$(CC) $(CFLAGS) -o $(TARGET) main.o 

clean: main
	$(RM) main.o point_mass_1D.o
run:
	$(RUN)$(TARGET)

plots: main run 
	$(R) $(RSCRIPT)
