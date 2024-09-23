import os

# Energy of the ion in keV
e = 20

for i in range( 100 ):
    os.system(f'./run.sh {i} {e}')