# radiometric_measurement

---

The module` radiometric_measurement `contains the class `ScanSurface`, which can be run in an ITOM environment to perform measurements with a radiometric scanning device developed at the Institute of Chemical Engineering, Ulm University.
    
The class ScanSurface controlling the measurement device is connected to the spectrometer, the power supply and the tinkerforge bricks, which control the stepper motors and end-stop switches. Must be run in an ITOM python environment.

## Usage
The following cole block shows an example of a script initializing an instance of the class `ScanSurface` class and starting a measurement. This script must be executed in ITOM. In line 3, the class instance `sc` gets initialized with the parameters described in the following section. After connecting to the hardware in line 8, the measurement is started in line 10 after a waiting time to leave the lab (line 9).

	from Skripts.radiometric_measurement import ScanSurface
	import time
	
	sc = ScanSurface(reactor='example_name', comment='example_comment', measurements=5, power_supply=True, current=0.5,
	                 slicelist=[10.8, 11.19, 11.58, 11.97, 12.36], 
	                 int_time_list=[0.142, 0.152, 0.163, 0.174, 0.185]
	                 )  # initialize an instance of the ScanSurface class
	sc.connect()
	time.sleep(2)  # an apropirate time to leave the lab
	sc.start_measurement()  # start the measurement

### Attributes constructed by \_\_init\_\_:

**reactor** ***(str):*** The name of the experiment, will be used for output folder and files.

**comment** ***(str):*** A comment for additional information, written to the logfile.

**measurements** ***(int):*** The size of measurements of a quadratic canvas in x/y-direction. The points are evenly spaced by 3.9 mm. The maximum total size is 153 (59.67 cm).

**power\_supply** ***(bool):*** If `True`, the script controlls the power supply with the given parameters for voltage and current. If `False`, a dark spectrum for the chosen integration time must be provided in the format "dark_'integration time'\_av'averages'.csv", whereas 'integration time' and 'averages' must be replaced with the respective numbers. The provided python script *darkspec.py* can help with that.

**voltage** ***(number, optional, default = 5.0):*** The voltage set for the power supply in Volt.

**current** ***(number, optional, default = 0.7):*** The current set for the power supply in Ampere.

**slicelist** ***(list):*** Sets up a list for all z-positions, in cm, to be measured: `[pos1, pos2, pos3...]`.

**int\_time\_list** ***(list):*** Sets up a list of the integration time for all z-positions, in seconds: `[time1, time2, time3...]`. Must be of the same length as parameter slicelist.

**averages** ***(int, optional, default = 2):*** The number of spectrometer measurements which result in one averaged measurement sent to python. A new transfer curve is needed for different average settings.

**x_start** ***(number, optional, default = 0.0):*** The x-position of the first measured pixel in cm.

**y_start** ***(number, optional, default = 0.0):*** The y-position of the first measured pixel in cm.

**log\_output** ***(boolean, optional, defalut = True):*** Saves a logfile, this is mandatory for further evaluation.

**linspace** ***(tuple, optional, default = False):*** Sets up a slicelist in the manner of the numpy linspace function *(start, stop, num)*. Linspace definitions can't be defined together with a list definition (slicelist or int\_time\_list)!

**int\_time\_linspace** ***(tuple, optional, default = False):*** Sets up an integration time list in the manner of the numpy linspace function *(start, stop, num)*. Must be used together with parameter linspace.

**sendmail** ***(bool, default is False):*** If `True`, the script will send a status mail after the 50th measurement of a respective z-slice. Attributes should be changed in the code.


### Attributes for status and debugging:

**is\_homed** ***(bool):*** Indicates overall homing status.

**bool\_home\_x** ***(bool):*** Indicates x homing status.

**bool\_home\_y** ***(bool):*** Indicates y homing status.

**z\_is\_homed** ***(bool):*** Indicates z homing status.

**status** ***(str):*** Indicated status `'idle'`, `'measuring'` or `'homing'`.

**path** ***(str):*** The filepath for result export. Normally consisting of '02Results/' the time as "%y%m%d_%H%M%S" and the string of the attribute reactor.

#### Attributes for the status mail function:

If sendmail is `True`, the script will send a status mail after the 50th measurement of a respective z-slice. Attributes should be changed in the code ([line 346](./Skripts/radiometric_measurement.py#L347-348)):

**yag** = yagmail.SMTP('fill your email here', 'fill your password here or use the keychain functionality of yagmail')

**mailto** ***(str):*** The receivers e-mail address.

**mailsubj** ***(str):*** Automatically generated subject.

**mailcont** ***(str):*** Automatically generated content.

**approx\_string** ***(str):*** Automatically generated string of the approximate ending time of the scan.


#### Attributes for Tinkerforge and Steppermotor setup:

For setting up a new machine, or if tinkerforge hardware is replaced, the respetive UIDs must be typed in the code ([line 300](./Skripts/radiometric_measurement.py#L300-304)):
	
**UIDSTEP1** ***(str):*** UID of Stepper Brick 1
	
**UIDSTEP2** ***(str):*** UID of Stepper Brick 2 Since a core xy-setup is used, there is no such a x or y stepper.
	
**UIDSTEP3** ***(str):*** UID of Stepper Brick 3 (z-direction).
	
**UIDIO** ***(str):*** UID of IO-16 Bricklet.
	
**UIDIDO4** ***(str):*** UID of IO-4 Bricklet.

**Optional tinkerforge attributes:**

**HOST** ***(str):*** "localhost"
	
**PORT** ***(int):*** 4223
	
**max\_vel** ***(int):*** Velocity of stepper 1 and 2.
	
**ramping** ***(int):*** Ramping of stepper 1 und 2.
	
**z\_max\_vel** ***(int):*** Velocity of z-steppers.
	
**z\_ramping** ***(int):*** Ramping of z-steppers.


#### Methods:

Please refer to the detailed docstrings insde the module's [code](./Skripts/radiometric_measurement.py).

- connect()

- start\_measurement()

- disconnect()

##### Debugging methods for idle mode: 
(idle mode = connected and homed, no measurement running, `mode=='idle'`)

- get\_position()

- get\_z\_position()

- goto\_position(position)

- tune\_inttime()

- manual\_homing()

- reset\_queue()


