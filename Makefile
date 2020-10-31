NAME = main

CC = g++
FLAGS = -Wall -Werror -O3
LIBRARIES = -lm
INCLUDES = -I.

HEADERS = analytic.h
SOURCES = main.cpp analytic.cpp

OBJECTS_DIRECTORY = objects/
OBJECTS_LIST = $(SOURCES:.cpp=.o)
OBJECTS = $(addprefix $(OBJECTS_DIRECTORY), $(OBJECTS_LIST))

GREEN = \033[0;32m
RED = \033[0;31m
RESET = \033[0m

.PHONY: all clean fclean re

all: $(NAME)


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

	