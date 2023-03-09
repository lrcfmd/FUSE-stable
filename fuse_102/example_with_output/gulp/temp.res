# 
# Keywords:
# 
opti conp c6  
# 
# Options:
# 
title
ASE calculation                                                                 
end
cell
   3.955151   3.955121   3.955194  90.000021  90.000010  90.000440
fractional  1   
Sr    core 0.0053686 0.9985056 0.2570004 2.00000000 1.00000 0.00000             
O     core 0.5053703 0.9984995 0.7570005 -2.0000000 1.00000 0.00000             
O     core 0.5053705 0.4985007 0.2569981 -2.0000000 1.00000 0.00000             
Ti    core 0.5053712 0.4985007 0.7569983 4.00000000 1.00000 0.00000             
O     core 0.0053715 0.4984971 0.7569992 -2.0000000 1.00000 0.00000             
totalenergy          -158.5835462535 eV
species   3
Sr     core    2.000000                  
O      core   -2.000000                  
Ti     core    4.000000                  
buck     
O     core Sr    core  1952.39000     0.336850  19.22000      1.20 12.00
buck     
O     core Ti    core  4590.72790     0.261000  0.000000      1.20 12.00
buck     
O     core O     core  1388.77000     0.362620  175.0000      1.20 12.00
lbfgs_order 5000
time        300.0
maxcyc opt     1500
maxcyc fit     1500
dump temp.res                                                    
