NAME = main

CC = g++
polus: CC = -mpixlC
bg: CC = mpixlcxx_r

# note that '-fopenmp' is a flag and not a library
FLAGS = -Wall -O0 -g -fopenmp -D FAKEMPI
polus: FLAGS = -qsmp=omp
bg: FLAGS = -qsmp=omp


LIBRARIES = -lm
INCLUDES = -I.

HEADERS = problem.h comm.hpp block.hpp
SOURCES = main.cpp comm.cpp block.cpp 

OBJECTS_DIRECTORY = objects/
OBJECTS_LIST = $(SOURCES:.cpp=.o)
OBJECTS = $(addprefix $(OBJECTS_DIRECTORY), $(OBJECTS_LIST))

GREEN = \033[0;32m
RED = \033[0;31m
RESET = \033[0m

.PHONY: all clean fclean re

all: $(NAME)
polus: modules $(NAME)
bg: $(NAME)

modules:
	module load SpectrumMPI

$(NAME): $(OBJECTS_DIRECTORY) $(OBJECTS)
	@$(CC) $(FLAGS) $(LIBRARIES) $(INCLUDES) $(OBJECTS) -o $(NAME)
	@echo -e "$(GREEN)$(NAME) was created$(RESET)"

$(OBJECTS_DIRECTORY):
	@mkdir -p $(OBJECTS_DIRECTORY)
	@echo -e "$(GREEN)$(OBJECTS_DIRECTORY) was created$(RESET)"

$(OBJECTS_DIRECTORY)%.o: %.cpp $(HEADERS)
	@$(CC) $(FLAGS) -c $(INCLUDES) $< -o $@
	@echo -e "$(GREEN)\t" $@ " was created$(RESET)"

clean:
	@rm -rf $(OBJECTS_DIRECTORY)
	@echo -e "$(NAME): $(RED)$(OBJECTS_DIRECTORY) was deleted$(RESET)"

fclean: clean
	@rm -f $(NAME)
	@echo -e "$(NAME): $(RED)$(NAME) was deleted$(RESET)"

re:
	@$(MAKE) fclean
	@$(MAKE) all

	