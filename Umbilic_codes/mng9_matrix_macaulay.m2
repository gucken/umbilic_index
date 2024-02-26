-- Macaulay2 script to compute primary decomposition of the determinant of the matrix M below
-- Execute with Macaulay command input "mng9_matrix_macaulay.m2"
-- John Guckenheimer, January 5, 2024

R = QQ[y1,y2,y3]

M = matrix{
{56, 42*y1, - 2 + 30*y1^2, - 6*y1 + 20*y1^3, - 12*y1^2 + 12*y1^4, - 20*y1^3 + 6*y1^5, - 30*y1^4 + 2*y1^6, -42*y1^5, -56*y1^6},
{56, 42*y2, - 2 + 30*y2^2, - 6*y2 + 20*y2^3, - 12*y2^2 + 12*y2^4, - 20*y2^3 + 6*y2^5, - 30*y2^4 + 2*y2^6, -42*y2^5, -56*y2^6},
{56, 42*y3, - 2 + 30*y3^2, - 6*y3 + 20*y3^3, - 12*y3^2 + 12*y3^4, - 20*y3^3 + 6*y3^5, - 30*y3^4 + 2*y3^6, -42*y3^5, -56*y3^6},
{56, 0, -2 , 0, 0, 0, 0, 0, 0},
{ 0,    -7,        -12*y1,         -15*y1^2,            -16*y1^3,           -15*y1^4,           -12*y1^5,  -7*y1^6,        0},
{ 0,    -7,        -12*y2,         -15*y2^2,            -16*y2^3,           -15*y2^4,           -12*y2^5,  -7*y2^6,        0},
{ 0,    -7,        -12*y3,         -15*y3^2,            -16*y3^3,           -15*y3^4,           -12*y3^5,  -7*y3^6,        0},
{ 0,    -7,        0,  0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 1, 1, 1}}

Md = det M;

mdd = degree Md

Mdt = terms Md;

Mdtn = #Mdt

Mdi = ideal Md;

Mdip = primaryDecomposition Mdi;

Mdipn = #Mdip

Mdip_{0..5}

Mdip6 = Mdip_6;

pr = Mdip6_0;

prt = terms pr;

prtn = #prt

prd = degree pr
