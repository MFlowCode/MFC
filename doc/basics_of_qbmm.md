'qbmm'               : 'T',
    * Turns on QBMM. It uses CHyQMOM.
'nnode'              : 4,
    * Number of quadrature points. 
    * Currently 4 is the only option.
'dist_type'          : 2,
    * 1 = Binormal distribution
    * 2 = Lognormal in R, Normal in V (Rdot)
    * Both of these use mean in R = 1, mean in V = 0
'sigR'               : 0.1,
    * Standard deviation of distribution in R
'sigV'               : 0.1,
    * Standard deviation of distribution in V (Rdot)
'rhoRV'              : 0.0,
    * Initial covariance between R and V directions (usually 0)
