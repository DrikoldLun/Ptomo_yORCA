import obspy
import glob, os

path = 'Pdata/vertical'
ch = 'Z'
filenames = glob.glob(os.path.join(path,'*/*.Z'))
for filename in filenames:
    st = obspy.read(filename)[0]
    if st.stats.channel != ch:
        st.stats.channel = ch
        st.write(filename,'SAC')
        print(filename)
    else:
        continue
