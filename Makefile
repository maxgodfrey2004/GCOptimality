.PHONY: data run all

REXE := Rscript.exe
CC := g++

DATA_FILE_NAME := aadata.txt

all: similarity.out 

run: similarity.out
	./similarity.out $(DATA_FILE_NAME)

similarity.out: random_similarity.cpp
	$(CC) random_similarity.cpp -o similarity.out

data:
	$(REXE) gen_aadata.r $(DATA_FILE_NAME)
