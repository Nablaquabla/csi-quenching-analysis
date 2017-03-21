CsI quenching factor measurements at TUNL
c =============================================================================
c                           Cell Cards
c =============================================================================
 1  0  99  $ outer world
c ====== CsI Detector =====
 2  2  -4.5100 (-2 3 -1) $ CsI crystal
 3  3  -0.8300 (-2 3 1 -4) : (-3 5 -4) $ Mylar coating
 4  4  -2.6989 (-2 5 4 -6):(-5 7 -6):(-10 6 -2 8):(-10 6 -9 7)$ Al housing  
c ====== Eljen 299 PSD Detector ======
 5  5  -1.08 (-21 22 -23) $ EJ 299-33A
 6  4  -2.6989 (-25 24 -23) #5 $ Al housing
c ====== Aluminum Table =====
 7  4  -2.6989 (-30 32 -31)
c ====== Air ======
 8  1  -0.00125 -98 
 90 1  -0.00125 -99 #2 #3 #4 #5 #6 #7 #8
 
c =============================================================================
c                           Surface Cards
c =============================================================================
c --- CsI Target ---
1 cz 0.95
2 pz 2
3 pz -3.1
4 cz 1.04
5 pz -3.2
6 cz 1.12
7 pz -3.3
8 pz 1.6
9 pz -2.3
10 cz 1.27
c --- Plastic Scintillator ---
21 1 cx 5.715
22 1 px 0
23 1 px 4.76
24 1 px -0.25
25 1 cx 5.9615
c --- Aluminum Table ---
30 cz 85.73
31 pz -42.13
32 pz -43.4
c --- Cookie Cutter for SDEF ---
98 rec -5.5 0 -0.65  1 0 0  0 2.7 0  0 0 4.4
c --- World ---
99 so 150.0

c =============================================================================
c                           Data Cards
c =============================================================================
imp:n 0 1 1 1 1 1 1 1 1
mode n
c ======================== Rotation Cards =====================================    
tr1 75.5231 48.9516 0 0.839146 0.543907 0 -0.543907 0.839146 0 0 0 1
c ================================== Materials ================================
c ----------- Air ---------
m1   8016.70c 23 &
     7014.70c 75 &
     18040.70c 1 $
c --- CsI @ 4.51 g/cc ---
m2   55133.70c 1 &
     53127.70c 1 $
c --- Mylar Reflector @ 0.83 g/cc ---
m3    6000.70c 11 &
      8016.70c 5 &
      1001.70c 12 $
c --- Aluminum 6061 Alloy @ 2.6989 g/cc ---
m4   13027.70c 96.55 &
     25055.70c 0.15 &
     29063.70c 0.21 &
     29065.70c 0.09 &
     30000.70c 0.25 &
     24050.70c 0.01 &
     24052.70c 0.16 &
     24053.70c 0.02 &
     24054.70c 0.01 &
     26054.70c 0.04 &
     26056.70c 0.64 &
     26057.70c 0.02 &
     40090.70c 0.05 &
     40091.70c 0.01 &
     40092.70c 0.02 &
     40094.70c 0.02 &
     22046.70c 0.01 &
     22047.70c 0.01 &
     22048.70c 0.11 &
     22049.70c 0.01 &
     22050.70c 0.01 &
     12024.70c 0.8 &
     12025.70c 0.1 &
     12026.70c 0.1 &
     14028.70c 0.55 &
     14029.70c 0.03 &
     14030.70c 0.02
c --- Plastic Scintillator EJ-299-33A @ 1.08 g/cc ---
m5    6000.70c 0.514 &
      1001.70c 0.486 $
c ================================== Source ===================================
c  Monochromatic Neutron Source @ XYZ MeV
c =============================================================================
SDEF POS = -5 0 0 &
     AXS = 1 0 0 &
     VEC = 1 0 0 &
     DIR = 1 &
     PAR = 1 &
     RAD = d1 &
     ERG = d2 &
     CEL = 8
SI1 0 4.4
SP1 -21 1
SI2 a 3.071 3.105 3.139 3.172 3.205 3.236 3.267 3.298 3.328 3.357 3.386 &
      3.415 3.443 3.471 3.499 3.526 3.554 3.580 3.607 3.634 3.660 3.686 & 
      3.712 3.737 3.763 3.788 3.813 3.838 3.863 3.888 3.912 3.937 3.961 & 
      3.985 4.009 4.033 4.057 4.081 4.104 4.128 4.151 4.175 4.198 4.221 &
      4.244 4.267 4.290 4.313 4.336 4.359
SP2   0 0 48 603 2024 3907 5870 7634 9003 9707 10082 &
      10893 12426 14041 15842 18115 21089 23793 26500 29817 33474 37386 &
      41187 45185 49024 52500 55113 56040 55230 51818 45933 39136 31674 & 
      23424 15011 7918 2941 252 0 0 0 0 0 0 &
      0 0 0 0 0 0
c ================================== Tallies ==================================
F1:n 22
c ================================== Misc =====================================
IPOL 0 0 0 0 J 2 2 2 5
FILES 21 DUMN1
NPS   1.0E9
