import os

ViscousJobs = [ ['C79E1D3C', 'CD9D3050'], \
                ['9CB03CEF', 'D6BAC936', 'FB822062', 'B3AAC9C8'], \
                ['1C0780C8', '2060F55A', '07C33719', '939D6718']]

MonopoleJobs = [['EF54219C', '7FC6826B', '59D05DE9', '6B22A317', 'AF46C382', '46AA7AF8'], \
                ['34DBFE14', '14B6198D', 'AE37D842', 'CC4F7C44', '66CFF8CC', '303B925A'], \
                ['AD63A4A5', '8A341282', '622DEC78', 'AB04C64D', '63850240' ]]


itr = 10
vTimes = [[0 for x in range(10)] for y in range(itr)]
mTimes = [[0 for x in range(17)] for y in range(itr)]  

# Viscous Tests
for k in range(0,itr):
    os.system("rm -rf RunTimeTests/ViscousTests/*")
    for i in ViscousJobs:
        for j in i:
            os.system("cp tests/" + j + "/case.py RunTimeTests/ViscousTests/case" + j + ".py")
            os.system("./mfc.sh run RunTimeTests/ViscousTests/case" + j + ".py")
    
    f = open("RunTimeTests/ViscousTests/time_data.dat")
    print(f)
    for m in range(0, 10):
        vTimes[k][m] = float(f.readline().split()[1])

meanVTimes = []
for m in range(0,10):
    testSum = 0
    for j in range(0, itr):
        testSum = testSum + vTimes[j][m]
    meanVTimes.append(testSum/itr)

# Monopole Tests
for k in range(0,itr):
    os.system("rm -rf RunTimeTests/MonopoleTests/*")
    for i in MonopoleJobs:
        for j in i:
            os.system("cp tests/" + j + "/case.py RunTimeTests/MonopoleTests/case" + j + ".py")
            os.system("./mfc.sh run RunTimeTests/MonopoleTests/case" + j + ".py")
    
    f = open("RunTimeTests/MonopoleTests/time_data.dat")
    for m in range(0, 17):
        mTimes[k][m] = float(f.readline().split()[1])   

meanMTimes = []
for m in range(0,17):
    testSum = 0
    for j in range(0, itr):
        testSum = testSum + mTimes[j][m]    
    meanMTimes.append(testSum/itr)
        
f = open("RunTimeTests/MeanTimeData", 'x')
f.write("Viscous\n")
for i in range(0,10):
    f.write(str(meanVTimes[i]) + "\n")
f.write("Monopole\n")
for i in range(0,17):
    f.write(str(meanMTimes[i]) + "\n")
f.close()