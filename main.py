import deep


alpha=30
border=20
a=2.46
mixshift=[0.1,0.1,0]
thickness=3.17
actual_part=1
excication_part=0.00025
bounding_border=0.75


ROOT_DIR = os.path.dirname(os.path.abspath(__file__)) 
sys.path.append(ROOT_DIR)


if os.path.exists("storage"): 
    shutil.rmtree("storage")
    # Создаем папку
os.makedirs("storage")

   

for j in range(7):

    mixshift=np.array([np.random.uniform(low=-a/math.sqrt(3), high=a/math.sqrt(3), size=2)[0], np.random.uniform(low=-a/math.sqrt(3), high=a/math.sqrt(3), size=2)[1], 0]) 
    

    subname_m=ROOT_DIR+"/storage/logs_m_"+str(j)
    subname_s=ROOT_DIR+"/storage/logs_s_"+str(j)


    x=bigraphene(alpha, border, a, mixshift, thickness, bounding_border)
    x.ShowPlot()

    x.StartModelling(bounding_border, ROOT_DIR, ROOT_DIR+"/opt.xyz")
    
    x.ReadEnergy(ROOT_DIR+"/out")
    
    x.MoveToCSVLogs(subname_m+".csv", 0)
    x.DumpToXYZ(subname_m+".xyz")


    
    for i in range(300):
        
        
        x.AddHydrogenPartly(excication_part, bounding_border)   
        x.StartModelling(bounding_border, ROOT_DIR, ROOT_DIR+"/opt.xyz")
        x.ReadEnergy(ROOT_DIR+"/out")

        n=np.shape(x.hydrogen_low)[0]+np.shape(x.hydrogen_high)[0]
        x.MoveToCSVLogs(subname_m+".csv", n)

        x.DumpToXYZ(subname_m+".xyz")

        print("Superteration "+str(i)+" done.")
        
        # if i//10 == 0:
        #     os.makedirs("storage/serg_"+str(i))
        #     dest_directory =ROOT_DIR+ "/storage/serg_"+str(i)


        #     file_name = "out"
        #     copy_files_with_name(ROOT_DIR, dest_directory, file_name)
        #     file_name = "opt.xyz"
        #     copy_files_with_name(ROOT_DIR, dest_directory, file_name)
