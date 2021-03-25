NAME = task

HEADERS = problem.h comm.hpp block.hpp
SOURCES = main.cpp comm.cpp block.cpp 
CU_SOURCES = mult.cu


# bsub -n 20 -oo cuda20.log -gpu "num=1:mode=shared:j_exclusive=yes" mpiexec ./task
CC = mpixlC
CU_CC = nvcc
# polus: CC = -mpixlC
# bg: CC = mpixlcxx_r
# bg_: CC = mpixlcxx_r

# mpisubmit.pl -p 10 --stdout task-256-10.log task -- -b 1 2 5 -n 256 -k 32 -t 0.0625 -l 3.1415 3.1415 3.1415
# note that '-fopenmp' is a flag and not a library
FLAGS = -Wall -D FAKEMPI -qsmp=omp
polus: FLAGS = -qsmp=omp
# bg: FLAGS = -qsmp=omp
# # No opemMP version
# bg_: FLAGS = 

LIBRARIES = -lm -L/usr/local/cuda/lib64 -lcudart
INCLUDES = -I.

OBJECTS_DIRECTORY = objects/cpp/
CUBJECTS_DIRECTORY = objects/cu/

OBJECTS_LIST = $(SOURCES:.cpp=.o)
CUBJECTS_LIST = $(CU_SOURCES:.cu=.o)

OBJECTS = $(addprefix $(OBJECTS_DIRECTORY), $(OBJECTS_LIST))
CUBJECTS = $(addprefix $(CUBJECTS_DIRECTORY), $(CUBJECTS_LIST))

GREEN = \033[0;32m
RED = \033[0;31m
RESET = \033[0m

.PHONY: all clean fclean re

all: $(NAME)
polus: $(NAME)
# bg: $(NAME)

# bg_: $(NAME)

# modules:
# 	module load SpectrumMPI

$(NAME): $(OBJECTS_DIRECTORY) $(CUBJECTS_DIRECTORY) $(OBJECTS) $(CUBJECTS)
	@$(CC) $(FLAGS) $(LIBRARIES) $(INCLUDES) $(OBJECTS) $(CUBJECTS) -o $(NAME)
	@echo -e "$(GREEN)$(NAME) was created$(RESET)"

$(OBJECTS_DIRECTORY):
	@mkdir -p $(OBJECTS_DIRECTORY)
	@echo -e "$(GREEN)$(OBJECTS_DIRECTORY) was created$(RESET)"

$(CUBJECTS_DIRECTORY):
	@mkdir -p $(CUBJECTS_DIRECTORY)
	@echo -e "$(GREEN)$(CUBJECTS_DIRECTORY) was created$(RESET)"

$(OBJECTS_DIRECTORY)%.o: %.cpp $(HEADERS)
	@$(CC) $(FLAGS) -c $(INCLUDES) $< -o $@
	@echo -e "$(GREEN)\t" $@ " was created$(RESET)"

$(CUBJECTS_DIRECTORY)%.o: %.cu $(HEADERS)
	@$(CU_CC) -c $(INCLUDES) $< -o $@
	@echo -e "$(GREEN)\t" $@ " was created$(RESET)"

clean:
	@rm -rf $(OBJECTS_DIRECTORY)
	@echo -e "$(NAME): $(RED)$(OBJECTS_DIRECTORY) was deleted$(RESET)"
	@rm -rf $(CUBJECTS_DIRECTORY)
	@echo -e "$(NAME): $(RED)$(CUBJECTS_DIRECTORY) was deleted$(RESET)"

fclean: clean
	@rm -f $(NAME)
	@echo -e "$(NAME): $(RED)$(NAME) was deleted$(RESET)"

re:
	@$(MAKE) fclean
	@$(MAKE) all

