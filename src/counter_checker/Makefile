.PHONY: all clean

all: counter_checker_rf41 counter_checker_rf32

counter_checker_rf41: counter_checker.c
	gcc -o counter_checker_rf41 -O3 -Wall -I ../../include -D RF41 -D STANDALONE_SETUP counter_checker.c

counter_checker_rf32: counter_checker.c
	gcc -o counter_checker_rf32 -O3 -Wall -I ../../include -D RF32 -D STANDALONE_SETUP counter_checker.c

clean:
	rm counter_checker_rf32 counter_checker_rf41
