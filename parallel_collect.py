import os
if not os.path.exists('log'):
    os.mkdir('log')
#mseed_path = '/home/lun/scratch-lun/OBSdata/yORCA_data/Mseed'
stalst = os.listdir('../../../OBSdata/yORCA_data/Mseed')
for sta in stalst:
    run_log = 'log/run_'+sta+'.log'
    run_err = 'log/run_'+sta+'.err'
    os.system('nohup python -u collect_P.py -S'+sta+' -Call paraP.cfg 1> '+run_log+' 2> '+run_err+' &')
