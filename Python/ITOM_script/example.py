from Skripts.radiometric_measurement import ScanSurface
import time

sc = ScanSurface(reactor='example_name', comment='example_comment', measurements=5, power_supply=True, current=0.5,
                 slicelist=[10.8, 11.19, 11.58, 11.97, 12.36], 
                 int_time_list=[0.142, 0.152, 0.163, 0.174, 0.185]
                 )  # initialize an instance of the ScanSurface class
sc.connect()
time.sleep(2)  # an apropirate time to leave the lab
sc.start_measurement()  # start the measurement