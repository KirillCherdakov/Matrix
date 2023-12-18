.PHONY: all s21_matrix.a test gcov_report style_check memory_check valgrind clean_gcov clean

CC=gcc
C_FLAGS=-Wall -Wextra -Werror -std=c11 -pedantic
MAC_FLAGS=-lcheck -lm
LINUX_FLAGS=-lcheck -lsubunit -lm -lrt -lpthread
GCOV_FLAGS=-fprofile-arcs -ftest-coverage
VALGRIND_FLAGS_v2 = --log-file="valgrind.txt" --track-origins=yes --trace-children=yes --leak-check=full --leak-resolution=med --track-fds=yes --show-leak-kinds=all
VALGRIND_FLAGS=--leak-check=full

FILES=$(wildcard *.c)
OBJ=$(patsubst %o, %o, $(patsubst %.c,%.o,$(FILES)))
TEST_FILES = $(wildcard unit_tests/s21_matrix_tests.c)

OS := $(shell uname -s)

ifeq ($(OS),Darwin)
	OPEN=open
endif
ifeq ($(OS),Linux)
	OPEN=xdg-open
endif

all: s21_matrix.a gcov_report

s21_matrix.a: $(FILES)
	$(CC) -c $(C_FLAGS) $(FILES)
	ar rc $@ $(OBJ)
	ranlib $@
	rm -rf *.o

test: s21_matrix.a
ifeq ($(OS),Darwin)
	$(CC) $(C_FLAGS) $(TEST_FILES) $< -o test $(MAC_FLAGS)
else
	$(CC) $(C_FLAGS) $(TEST_FILES) $< -o test $(LINUX_FLAGS)
endif
	./test

gcov_report: clean
ifeq ($(OS),Darwin)
	$(CC) $(C_FLAGS) $(GCOV_FLAGS) $(FILES) $(TEST_FILES) -o test $(MAC_FLAGS)
else
	$(CC) $(C_FLAGS) $(GCOV_FLAGS) $(FILES) $(TEST_FILES) -o test $(LINUX_FLAGS)
endif
	./test
	lcov -t "s21_test" -o postreport.info -c -d .
	genhtml -o report postreport.info
	$(OPEN) report/index.html

style_check:
	clang-format -i -style=google unit_tests/*.c unit_tests/*.h *.c *.h
	clang-format -n -style=google unit_tests/*.c unit_tests/*.h *.c *.h

memory_check: s21_matrix.a
	$(CC) $(C_FLAGS) $(TEST_FILES) $< -o memory_test $(LINUX_FLAGS)
	CK_FORK=no valgrind $(VALGRIND_FLAGS) ./memory_test

valgrind:
	valgrind $(VALGRIND_FLAGS_v2) ./test
	grep -n "ERROR SUMMARY" valgrind.txt | grep -v "0 errors"

clean_gcov:
	rm -rf *.gcda
	rm -rf *.gcno
	rm -rf test
	rm -rf report
	rm -rf postreport.info

clean: clean_gcov
	rm -rf *.o
	rm -rf *.a
	rm -rf test
	rm -rf memory_test