import ctout_read as ctor

entries=ctor.read('../test.ctout')

#print(entries)

for e in entries:
    print('Positive test time: ',e[0]/1440.,', time of symptoms appearance: ',e[1]/1440.,', ID: ',e[2],', Parent ID: ',e[3],', Number of successfully traced contacts: ',e[4])
