# Simple finite element assembler in J

This code assembles finite element matrices (related to bilinear forms) with P1 elements. The bilinear form can be an arbitrary combination of function values and their gradients. Writes output as triplets (value,i,j) into files (output.dat,outputi.dat,outputj.dat).
