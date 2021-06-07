import os,datetime
xmark = ['Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun']
#x mark
all_begin = datetime.date(2018, 4, 1)
k = 0 
linesx = []
for i in range(2018,2020):
    for j in range(1,13):
        if (i == 2018 and j < 4) or (i == 2019 and j > 6):
            continue
        hour_index = (datetime.date(i,j,1)-all_begin).days*24
        linesx.append(' '.join([str(hour_index),'a',xmark[k]])+'\n')
        k += 1
with open('x_mark.txt','w+') as f:
    f.writelines(linesx)
    f.close()

with open('query5.csv','r') as f:
    event_lst = f.readlines()[1:]
    f.close()

lines = []
for event in event_lst:
    event = event.split(',')
    mag = event[4]
    time = event[0].split('T')
    date = time[0].split('-')
    y, m, d = int(date[0]), int(date[1]), int(date[2])
    h = int(time[1].split(':')[0])
    minu = int(time[1].split(':')[1])
    hour_index = (datetime.date(y,m,d)-all_begin).days*24+h+minu/24.
    lines.append(' '.join([str(hour_index),mag,event[0],event[2],event[1],event[3]])+'\n')

outfile = 'event_h5.lst'
with open(outfile,'w+') as f:
    f.writelines(lines)
    f.close()

os.system('sort -n -k 1 '+outfile+' -o '+outfile)
