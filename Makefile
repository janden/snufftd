cc = gcc

libraries = -lm

executables = test_nufftd_spread \
              test_sub_snufftd_spread

test_nufftd_spread_sources = test_nufftd_spread.c nufftd_spread.c

test_sub_snufftd_spread_sources = test_sub_snufftd_spread.c sub_snufftd_spread.c

all: $(executables)

test_nufftd_spread: $(test_nufftd_spread_sources)
	$(cc) -o test_nufftd_spread $(test_nufftd_spread_sources) $(libraries)

test_sub_snufftd_spread: $(test_sub_snufftd_spread_sources)
	$(cc) -o test_sub_snufftd_spread $(test_sub_snufftd_spread_sources) $(libraries)

clean:
	rm -f $(executables)
