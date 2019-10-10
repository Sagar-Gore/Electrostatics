#!/usr/bin/python

## Electrostatics potential file parser ##
# Mapping electrostatic potentials at atomic positions (3D coordinates)
### usage: python potdx_parser_AQUASOL.py pdb_file potdx_file output_file

import sys

def potdx_parser(pdb_file, potdx_file, output_file):
    import sys
    import time
    start_time = time.time()
    pdb = open(pdb_file, "r")  # 3D atomic coordinates
    pot = open(potdx_file, "r")  # Electrostatic Potential file APBS format
    op = open(output_file, "w")  # Output File
    pdb_name = pdb_file[0:4]
    num = []
    potential = []
    grid_neighbours = {}
    count_line = 0
    at_cnt = 0
    count = 0
    grid_cnt = 0
    X_gp = []
    Y_gp = []
    Z_gp = []
    X_coord = []
    Y_coord = []
    Z_coord = []
    pyro_atom = []
    pyro_chain = []
    pyro_residue = []
    pyro_gnum = []
    atom_num = []
    atom_name = [] 
    residue = []
    chain = []
    res_seq_num = []

    # read pdb file and store 3D coordinates
    for line in pdb:
        if(line[0:6] == "ATOM  "):
            #print line.strip('\n')
            atom_num.append(int(line[6:11]))
            atom_name.append(line[12:16].strip())
            residue.append(line[17:26])
            chain.append(line[21:22])
            res_seq_num.append(int(line[22:26]))
            X_coord.append(float(line[30:38]))
            Y_coord.append(float(line[38:46]))
            Z_coord.append(float(line[46:54]))
            at_cnt = at_cnt + 1
 
    # read AQUASOL potdx file and store electrostatic potential values
    for line in pot:
        word = line.split()
        count_line = count_line + 1
        if (word[0] == 'object' and word[2] == 'class' and word[3] == 'gridpositions' and word[4] == 'counts'):
            X_points = int(word[5])
            Y_points = int(word[6])
            Z_points = int(word[7])
            count = X_points * Y_points * Z_points
        else:
            if (word[0] == 'origin'):
                X_Origin = float(word[1])
                Y_Origin = float(word[2])
                Z_Origin = float(word[3])
                #print "X_Origin=", X_Origin, "Y_Origin=", Y_Origin, "Z_Origin=", Z_Origin
            else:
                if (word[0] == 'delta' and (word[2] == '0.000000e+00' and word[3] == '0.000000e+00')):
                    X_GSpace = float(word[1])
                if (word[0] == 'delta' and (word[1] == '0.000000e+00' and word[3] == '0.000000e+00')):
                    Y_GSpace = float(word[2])
                if (word[0] == 'delta' and (word[1] == '0.000000e+00' and word[2] == '0.000000e+00')):
                    Z_GSpace = float(word[3])
                    #print "X_Grid=", X_GSpace, "Y_Grid=", Y_GSpace, "Z_Grid=", Z_GSpace
                if (count_line > 11 and count_line <= (count / 3) + 11):
                    x = line.strip('\n')
                    y = x.split()
                    xa = float(y[0])
                    ya = float(y[1])
                    za = float(y[2])
                    potential.extend([xa, ya, za])
                    
    # Define grid points and store grid neighbours
    for i in range(0, X_points):
        for j in range(0, Y_points):
            for k in range(0, Z_points):
                num.append(str( str(i) + ':' +  str(j) + ':' +  str(k)))
    
    soo_X = []
    soo_Y = []
    soo_Z = []

    # make a mesh with grid points
    for i in range(0, X_points):
        xp = float(X_Origin + (i * X_GSpace))
        X_gp.append(xp)
        soo_X.append(xp - X_Origin)
    for j in range(0, Y_points):
        yp = float(Y_Origin + (j * Y_GSpace))
        Y_gp.append(yp)
        soo_Y.append(yp - Y_Origin)
    for k in range(0, Z_points):
        zp = float(Z_Origin + (k * Z_GSpace))
        Z_gp.append(zp)
        soo_Z.append(zp - Z_Origin)

#look through 3D mesh and fit atoms in a mesh
    for cnt in range(0,at_cnt):
        x = str(int((X_coord[cnt] - X_Origin)/X_GSpace)) 
        y = str(int((Y_coord[cnt] - Y_Origin)/Y_GSpace))
        z = str(int((Z_coord[cnt] - Z_Origin)/Z_GSpace))
        pyro_atom.append(atom_name[cnt])
        pyro_residue.append(residue[cnt])
        pyro_chain.append(chain[cnt])
        pyro_gnum.append(str(x + ':' + y + ':' + z))
					
#get potentials at the nearest grid points of atoms
    for cnt in range(0,at_cnt):
        g_XYZ = []
        pot_sum = 0

        g_XYZ = pyro_gnum[cnt].split(':')
        g_XYZ[0]  = int(g_XYZ[0])
        g_XYZ[1]  = int(g_XYZ[1])
        g_XYZ[2]  = int(g_XYZ[2])

        grid_neighbours.update({pyro_gnum[cnt]:[[g_XYZ[0],g_XYZ[1],g_XYZ[2]],[g_XYZ[0]+1,g_XYZ[1],g_XYZ[2]],[g_XYZ[0],g_XYZ[1]+1,g_XYZ[2]],[g_XYZ[0],g_XYZ[1],g_XYZ[2]+1],[g_XYZ[0]+1,g_XYZ[1]+1,g_XYZ[2]],[g_XYZ[0]+1,g_XYZ[1],g_XYZ[2]+1],[g_XYZ[0],g_XYZ[1]+1,g_XYZ[2]+1],[g_XYZ[0]+1,g_XYZ[1]+1,g_XYZ[2]+1]]})			
        for i in range(len(grid_neighbours[pyro_gnum[cnt]])):
            item = grid_neighbours[pyro_gnum[cnt]][i]
            vol = str(item[0])+':'+str(item[1])+':'+str(item[2])
            index = num.index(vol)
            pot_sum = pot_sum + potential[index]
        pot_temp = round((pot_sum/8),2)
        op.write('\t'.join([pdb_name,pyro_chain[cnt],str(pot_temp),pyro_atom[cnt], pyro_residue[cnt]]))
        op.write("\n")
    op.close()
    print 'time', time.time() - start_time
potdx_parser(sys.argv[1], sys.argv[2], sys.argv[3])

