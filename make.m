if isoctave()
    mex nufftd_spread_mx.c nufftd_spread.c

    mex sub_snufftd_spread_mx.c sub_snufftd_spread.c
else
    mex -largeArrayDims nufftd_spread_mx.c nufftd_spread.c

    mex -largeArrayDims sub_snufftd_spread_mx.c sub_snufftd_spread.c
end

clear functions;
