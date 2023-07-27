import sys
from array import array
input=sys.argv[1]
SYMM=[]
atom=[]
x=array('f')
y=array('f')
z=array('f')
atomtype=array('i')
SFAC=[]
lines=open(input).readlines()
output=input.split('.') [0]+".cif"
#print input, output
#spgr=lines[0].split()[2].split('_') [1]
spgr=lines[0].split("_")[1].split()[0]
#spgr=lines[0].split()[2].split('_')[1]
len1=float(lines[1].split() [2])
len2=float(lines[1].split() [3])
len3=float(lines[1].split() [4])
ang1=float(lines[1].split() [5])
ang2=float(lines[1].split() [6])
ang3=float(lines[1].split() [7])
SYMM.append("X,Y,Z\n")
print "data_1"
print "_symmetry_space_group_name_H-M   '",spgr,"'"
 

for i in range(2,len(lines)):
	if "ZERR" in lines[i]:
		ZERR=int(lines[i].split() [1])
	elif "LATT" in lines[i]:
		LATT=int(lines[i].split() [1])
	elif "SYMM" in lines[i]:
#		SYMM.append(lines[i].split() [1])
		X= lines[i] [5:len(lines[i])].split(',') [0].strip()
		Y= lines[i].split(',') [1].strip()
		Z= lines[i].split(',') [2].strip()
		if "+" in X[0]:
			X=X[1:len(X)]
		if "+" in Y[0]:
			Y=Y[1:len(Y)]
		if "+" in Z[0]:
			Z=Z[1:len(Z)]
		symm="%s,%s,%s\n" % (X,Y,Z)
		symm=symm.replace(' ','')
		SYMM.append(symm)
#		SYMM.append(lines[i] [5:len(lines[i])])
	elif "SFAC" in lines[i]:
		for j in range(1,len(lines[i].split())):
			SFAC.append(lines[i].split() [j])
	elif "FVAR" in lines[i]:
		no_one_cares=1
	elif "FVAR" in lines[i]:
		no_one_cares=1
	elif "UNIT" in lines[i]:
		no_one_cares=1
	elif "END" in lines[i]:
		break
	else:
		atom.append(lines[i].split() [0])
		atomtype.append(int(lines[i].split() [1]))
		x.append(float(lines[i].split() [2]))
		y.append(float(lines[i].split() [3]))
		z.append(float(lines[i].split() [4]))
#		print atom[-1],atomtype[-1]
#if "P" in spgr:
if LATT>0:
	for j in range(len(SYMM)):
		X=SYMM[j].split(',') [0]
		Y=SYMM[j].split(',') [1]
		Z=SYMM[j].split(',') [2]
		if "-" in X:
			X=X[1:len(X)]
		else:
			X="-"+X
		if "-" in Y:
			Y=Y[1:len(Y)]
		else:
			Y="-"+Y
		if "-" in Z:
			Z=Z[1:len(Z)]
		else:
			Z="-"+Z
		SYMM.append("%s,%s,%s" % (X,Y,Z))
	if "C" in spgr[0]:
		for j in range(len(SYMM)):
			X=SYMM[j].split(',') [0]
			Y=SYMM[j].split(',') [1]
			Z=SYMM[j].split(',') [2]
			X=X+"+1/2" #might have to put in some logic for if there is already a 1/2 there but might not
			Y=Y+"+1/2"
			SYMM.append("%s,%s,%s" % (X,Y,Z))
if len(SYMM)==4:
	if "C" in spgr[0]:
		for j in range(len(SYMM)): #WHY IS THIS GETTING C2/c as well as 106? 
			X=SYMM[j].split(',') [0]
			Y=SYMM[j].split(',') [1]
			Z=SYMM[j].split(',') [2]
			X=X+"+1/2" #might have to put in some logic for if there is already a 1/2 there but might not
			Y=Y+"+1/2"
			SYMM.append("%s,%s,%s" % (X,Y,Z))
	
#		SYMM.append("-X,-Y,-Z\n")
print "loop_"
print "_symmetry_equiv_pos_site_id"
print "_symmetry_equiv_pos_as_xyz"
#print "1 x,y,z"
for i in range(len(SYMM)):
	print i+1,SYMM[i],
print "_cell_length_a                   ",len1
print "_cell_length_b                   ",len2
print "_cell_length_c                   ",len3
print "_cell_angle_alpha                ",ang1
print "_cell_angle_beta                 ",ang2
print "_cell_angle_gamma                ",ang3
print "loop_"
print "_atom_site_label"
print "_atom_site_type_symbol"
print "_atom_site_fract_x"
print "_atom_site_fract_y"
print "_atom_site_fract_z"
for i in range(len(atom)):
#	print atom[i],SFAC[atomtype[i]-1],x[i],y[i],z[i]
	print "%-4s %4s %10.6f %10.6f %10.6f" % (atom[i],SFAC[atomtype[i]-1],x[i],y[i],z[i])
print "#END"
 
#		SFAC[]=lines[i].split() [1:len(lines[i].split())]

#	print spgr


