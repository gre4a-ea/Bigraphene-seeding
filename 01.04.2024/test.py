
int_left=1
count=2000

file_content = """units           metal
newton          on
dimension       3
boundary        p p p
atom_style      atomic

processors * * 1

read_data coo

pair_style mlip mlip.ini 
pair_coeff * *

neigh_modify    every 1 delay 0 check yes
variable        e equal etotal

thermo 10
thermo_style custom step epair etotal cpu

group freeze id 1:"""+str(count)+"""
fix 4 freeze setforce 0.0 0.0 0.0 

dump            1 all xyz 50 opt.xyz
dump_modify     1 append no element C H

minimize        1.0e-10 1.0e-6 20000 200000
      
print "$e" """

# Открываем файл 'input' в режиме записи
with open('input', 'w') as f:
    # Записываем текст в файл
    f.write(file_content)


int* arr = new int[5];
unsigned element = 2;
arr += element